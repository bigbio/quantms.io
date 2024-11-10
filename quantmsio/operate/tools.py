import os
import re
from collections import defaultdict
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO
import ahocorasick
from quantmsio.core.common import FEATURE_SCHEMA, IBAQ_SCHEMA, IBAQ_USECOLS, PSM_SCHEMA
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.query import Query, map_spectrum_mz
from quantmsio.core.openms import OpenMSHandler
from quantmsio.utils.pride_utils import get_unanimous_name
from quantmsio.utils.file_utils import load_de_or_ae, read_large_parquet


def generate_features_of_spectrum(
    parquet_path: str,
    mzml_directory: str,
    output_folder: str,
    file_num: int,
    partitions: list = None,
):
    """
    parquet_path: parquet file path
    mzml_directory: mzml file directory path
    """
    pqwriters = {}
    pqwriter_no_part = None
    filename = os.path.basename(parquet_path)
    p = Query(parquet_path)
    for _, table in p.iter_file(file_num=file_num):
        refs = table["reference_file_name"].unique()
        mzml_handlers = {ref: OpenMSHandler() for ref in refs}
        table[["number_peaks", "mz_array", "intensity_array"]] = table[["reference_file_name", "scan"]].apply(
            lambda x: map_spectrum_mz(
                x["reference_file_name"],
                x["scan"],
                mzml_handlers,
                mzml_directory,
            ),
            axis=1,
            result_type="expand",
        )
        if not partitions or len(partitions) > 0:
            for key, df in table.groupby(partitions):
                parquet_table = pa.Table.from_pandas(df, schema=PSM_SCHEMA)
                folder = [output_folder] + [str(col) for col in key]
                folder = os.path.join(*folder)
                if not os.path.exists(folder):
                    os.makedirs(folder, exist_ok=True)
                save_path = os.path.join(*[folder, filename])
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        else:
            parquet_table = pa.Table.from_pandas(table, schema=PSM_SCHEMA)
            if not os.path.exists(output_folder):
                os.makedirs(output_folder, exist_ok=True)
            save_path = os.path.join(*[output_folder, filename])
            if not pqwriter_no_part:
                pqwriter_no_part = pq.ParquetWriter(save_path, parquet_table.schema)
            pqwriter_no_part.write_table(parquet_table)
    if not partitions or len(partitions) == 0:
        if pqwriter_no_part:
            pqwriter_no_part.close()
    else:
        for pqwriter in pqwriters.values():
            pqwriter.close()


# gei unqnimous name
def map_protein_for_parquet(parquet_path, fasta, output_path, map_parameter, label):
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param parquet_path: psm_parquet_path or feature_parquet_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    :param label: feature or psm
    retrun: None
    """

    map_protein_names = defaultdict(set)
    if map_parameter == "map_protein_name":
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(name)
            map_protein_names[accession].add(name)
            map_protein_names[name].add(name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(accession)
            map_protein_names[accession].add(accession)
            map_protein_names[name].add(accession)
    change_and_save_parquet(parquet_path, map_protein_names, output_path, label)


def change_and_save_parquet(parquet_path, map_dict, output_path, label):
    pqwriter = None
    for df in read_large_parquet(parquet_path):
        df["protein_accessions"] = df["protein_accessions"].apply(lambda x: get_unanimous_name(x, map_dict))
        if label == "feature":
            schema = FEATURE_SCHEMA
        elif label == "psm":
            schema = PSM_SCHEMA
        parquet_table = pa.Table.from_pandas(df, schema=schema)
        if not pqwriter:
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()


def map_protein_for_tsv(path: str, fasta: str, output_path: str, map_parameter: str):
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param path: de_path or ae_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    map_protein_names = defaultdict(set)
    if map_parameter == "map_protein_name":
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(name)
            map_protein_names[accession].add(name)
            map_protein_names[name].add(name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(accession)
            map_protein_names[accession].add(accession)
            map_protein_names[name].add(accession)
    df, content = load_de_or_ae(path)
    df["protein"] = df["protein"].apply(lambda x: get_unanimous_name(x, map_protein_names))
    content += df.columns.str.cat(sep="\t") + "\n"
    for _, row in df.iterrows():
        content += "\t".join(map(str, row)).strip() + "\n"
    with open(output_path, "w", encoding="utf8") as f:
        f.write(content)


def get_modification_details(seq: str, mods_dict: dict, automaton: any, select_mods: list = None):
    if "(" not in seq:
        return (seq, [])
    total = 0
    modifications = []
    modification_details = []
    peptidoform = ""
    pre = 0
    for item in automaton.iter(seq):
        total += len(item[1])
        modification = item[1][1:-1]
        position = item[0] - total + 1
        name = re.search(r"([^ ]+)\s?", modification)
        name = name.group(1)
        if position == 0:
            peptidoform = f"[{name}]-{peptidoform}"
        elif item[0] + 1 == len(seq):
            peptidoform += seq[pre : item[0] - len(item[1]) + 1]
            peptidoform += f"-[{name}]"
        else:
            peptidoform += seq[pre : item[0] - len(item[1]) + 1]
            peptidoform += f"[{name}]"
        pre = item[0] + 1
        if modification in modifications:
            index = modifications.index(modification)
            modification_details[index]["fields"].append({"position": position, "localization_probability": 1.0})
        elif modification in select_mods:
            modifications.append(modification)
            modification_details.append(
                {
                    "modification_name": mods_dict[modification][0],
                    "fields": [{"position": position, "localization_probability": 1.0}],
                }
            )
    peptidoform += seq[pre:]
    return (peptidoform, modification_details)


def get_ahocorasick(mods_dict: dict):
    automaton = ahocorasick.Automaton()
    for key in mods_dict.keys():
        key = "(" + key + ")"
        automaton.add_word(key, key)
    automaton.make_automaton()
    return automaton


def get_field_schema(parquet_path):
    schema = pq.read_schema(parquet_path)
    return schema


PROTEIN_ACCESSION = r"\|([^|]*)\|"


def get_protein_accession(proteins: str = None):
    proteins = str(proteins)
    if "|" in proteins:
        return re.findall(PROTEIN_ACCESSION, proteins)
    else:
        return re.split(r"[;,]", proteins)


def transform_ibaq(df):
    def transform(row):
        map_dict = row["intensities"]
        return map_dict["sample_accession"], map_dict["channel"], map_dict["intensity"]

    df = df.explode("intensities")
    df.reset_index(drop=True, inplace=True)
    df[["sample_accession", "channel", "intensity"]] = df[["intensities"]].apply(
        transform, axis=1, result_type="expand"
    )
    df.drop(["intensities"], axis=1, inplace=True)
    return df


def genereate_ibaq_feature(sdrf_path, parquet_path):
    Sdrf = SDRFHandler(sdrf_path)
    sdrf = Sdrf.transform_sdrf()
    experiment_type = Sdrf.get_experiment_type_from_sdrf()
    p = Query(parquet_path)
    for _, df in p.iter_file(file_num=10, columns=IBAQ_USECOLS):
        df = transform_ibaq(df)
        if experiment_type != "LFQ":
            df = pd.merge(
                df,
                sdrf,
                left_on=["reference_file_name", "channel"],
                right_on=["reference", "label"],
                how="left",
            )
        else:
            df = pd.merge(
                df,
                sdrf,
                left_on=["reference_file_name"],
                right_on=["reference"],
                how="left",
            )
        df.drop(
            [
                "reference",
                "label",
            ],
            axis=1,
            inplace=True,
        )
        df["fraction"] = df["fraction"].astype(str)
        feature = pa.Table.from_pandas(df, schema=IBAQ_SCHEMA)
        yield feature


def write_ibaq_feature(sdrf_path, parquet_path, output_path):
    pqwriter = None
    for feature in genereate_ibaq_feature(sdrf_path, parquet_path):
        if not pqwriter:
            pqwriter = pq.ParquetWriter(output_path, feature.schema)
        pqwriter.write_table(feature)
    if pqwriter:
        pqwriter.close()
