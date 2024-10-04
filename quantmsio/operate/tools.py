import os
from collections import defaultdict
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO
from quantmsio.core.project import ProjectHandler
from quantmsio.core.feature import FEATURE_SCHEMA
from quantmsio.core.psm import PSM_SCHEMA
from quantmsio.operate.query import Query
from quantmsio.utils.pride_utils import get_unanimous_name
from quantmsio.utils.file_utils import read_large_parquet


def generate_features_of_spectrum(
    parquet_path: str,
    mzml_directory: str,
    output_path: str,
    label: str,
    file_num: int,
    partition: str = None,
):
    """
    parquet_path: parquet file path
    mzml_directory: mzml file directory path
    """
    pqwriters = {}
    pqwriter_no_part = None
    p = Query(parquet_path)
    for table in p.iter_file(file_num=file_num):
        table = p.inject_spectrum_msg(table, mzml_directory)
        if label == "feature":
            schema = FEATURE_SCHEMA
        else:
            schema = PSM_SCHEMA
        # save
        if partition == "charge":
            for key, df in table.groupby(["charge"]):
                parquet_table = pa.Table.from_pandas(df, schema=schema)
                save_path = output_path.split(".")[0] + "-" + str(key) + ".parquet"
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        elif partition == "reference_file_name":
            for key, df in table.groupby(["reference_file_name"]):
                parquet_table = pa.Table.from_pandas(df, schema=schema)
                save_path = output_path.split(".")[0] + "-" + str(key) + ".parquet"
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        else:
            parquet_table = pa.Table.from_pandas(table, schema=schema)
            if not pqwriter_no_part:
                pqwriter_no_part = pq.ParquetWriter(output_path, parquet_table.schema)
            pqwriter_no_part.write_table(parquet_table)
    # close f
    if not partition:
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
    from Bio import SeqIO

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
    for index, row in df.iterrows():
        content += "\t".join(map(str, row)).strip() + "\n"
    with open(output_path, "w", encoding="utf8") as f:
        f.write(content)


def load_de_or_ae(path):
    f = open(path, encoding="utf-8")
    line = f.readline()
    pos = 0
    content = ""
    while line.startswith("#"):
        pos = f.tell()
        content += line
        line = f.readline()
    f.seek(pos - 1)
    return pd.read_csv(f, sep="\t"), content


# register_file
def register_file_to_json(project_file, attach_file, category, replace_existing):
    register = ProjectHandler(project_json_file=project_file)
    register.add_quantms_file(attach_file, category, replace_existing)
    register.save_updated_project_info(output_file_name=project_file)
