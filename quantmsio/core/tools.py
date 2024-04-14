import base64
import json
import os
import string
from collections import defaultdict
from io import BytesIO

import ahocorasick
import matplotlib.pyplot as plt
import mygene
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import swifter
from Bio import SeqIO
from venn import venn

from quantmsio.core.feature import FeatureHandler
from quantmsio.core.openms import OpenMSHandler
from quantmsio.core.plots import plot_distribution_of_ibaq
from quantmsio.core.plots import plot_intensity_box_of_samples
from quantmsio.core.plots import plot_intensity_distribution_of_samples
from quantmsio.core.plots import plot_peptide_distribution_of_protein
from quantmsio.core.plots import plot_peptides_of_lfq_condition
from quantmsio.core.project import ProjectHandler
from quantmsio.core.psm import PSMHandler
from quantmsio.core.query import Parquet
from quantmsio.core.statistics import IbaqStatistics
from quantmsio.core.statistics import ParquetStatistics
from quantmsio.utils.pride_utils import generate_gene_name_map
from quantmsio.utils.pride_utils import get_gene_accessions
from quantmsio.utils.pride_utils import get_unanimous_name
from quantmsio.utils.report import report


# optional about spectrum
def map_spectrum_mz(mz_path: str, scan: str, mzml: OpenMSHandler, mzml_directory: str):
    """
    mz_path: mzML file path
    scan: scan number
    mzml: OpenMSHandler object
    """
    if mzml_directory.endswith("/"):
        mz_path = mzml_directory + mz_path + ".mzML"
    else:
        mz_path = mzml_directory + "/" + mz_path + ".mzML"
    mz_array, array_intensity = mzml.get_spectrum_from_scan(mz_path, int(scan))
    return mz_array, array_intensity, 0


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
    p = Parquet(parquet_path)
    for table in p.iter_file(file_num=file_num):
        table = p.inject_spectrum_msg(table, mzml_directory)
        if label == "feature":
            hander = FeatureHandler()
        else:
            hander = PSMHandler()
        # save
        if partition == "charge":
            for key, df in table.groupby(["charge"]):
                parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
                save_path = output_path.split(".")[0] + "-" + str(key) + ".parquet"
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        elif partition == "reference_file_name":
            for key, df in table.groupby(["reference_file_name"]):
                parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
                save_path = output_path.split(".")[0] + "-" + str(key) + ".parquet"
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        else:
            parquet_table = pa.Table.from_pandas(table, schema=hander.schema)
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


# option about gene


def generate_gene_acession_map(gene_names, species="human"):
    mg = mygene.MyGeneInfo()
    gene_accessions = mg.querymany(gene_names, scopes="symbol", species=species, fields="accession")
    gene_accessions_maps = defaultdict(list)
    for obj in gene_accessions:
        if "accession" in obj:
            gene_accessions_maps[obj["query"]].append(obj["accession"])
    return gene_accessions_maps


def map_gene_msgs_to_parquet(
    parquet_path: str,
    fasta_path: str,
    map_parameter: str,
    output_path: str,
    label: str,
    species: str,
):
    map_gene_names = generate_gene_name_map(fasta_path, map_parameter)
    pqwriter = None
    for df in read_large_parquet(parquet_path):
        df["gene_names"] = df["protein_accessions"].swifter.apply(lambda x: get_unanimous_name(x, map_gene_names))
        gene_names = list(set([item for sublist in df["gene_names"] for item in sublist]))
        gene_accessions_maps = generate_gene_acession_map(gene_names, species)
        df["gene_accessions"] = df["gene_names"].swifter.apply(lambda x: get_gene_accessions(x, gene_accessions_maps))
        if label == "feature":
            hander = FeatureHandler()
        elif label == "psm":
            hander = PSMHandler()
        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
            # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()


# plot venn
def plot_peptidoform_charge_venn(parquet_path_list, labels):
    data_map = {}
    for parquet_path, label in zip(parquet_path_list, labels):
        df = pd.read_parquet(parquet_path, columns=["peptidoform", "charge"])
        psm_message = "Total number of PSM for " + label + ": " + str(len(df))
        print(psm_message)
        unique_pep_forms = set((df["peptidoform"] + df["charge"].astype(str)).to_list())
        pep_form_message = "Total number of Peptidoform for " + label + ": " + str(len(unique_pep_forms))
        print(pep_form_message)
        data_map[label] = unique_pep_forms
    plt.figure(figsize=(16, 12), dpi=500)
    venn(
        data_map,
        legend_loc="upper right",
        figsize=(16, 12),
        fmt="{size}({percentage:.1f}%)",
    )
    plt.savefig("pep_form_compare_venn.png")


def plot_sequence_venn(parquet_path_list, labels):
    data_map = {}
    for parquet_path, label in zip(parquet_path_list, labels):
        sequence = pd.read_parquet(parquet_path, columns=["sequence"])
        unique_seqs = set(sequence["sequence"].to_list())
        pep_message = "Total number of peptide for " + label + ": " + str(len(unique_seqs))
        print(pep_message)
        data_map[label] = unique_seqs
    plt.figure(figsize=(16, 12), dpi=500)
    venn(
        data_map,
        legend_loc="upper right",
        figsize=(16, 12),
        fmt="{size}({percentage:.1f}%)",
    )
    plt.savefig("sequence_compare_venn.png")


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


def change_and_save_parquet(parquet_path, map_dict, output_path, label):
    pqwriter = None
    for df in read_large_parquet(parquet_path):
        df["protein_accessions"] = df["protein_accessions"].apply(lambda x: get_unanimous_name(x, map_dict))
        if label == "feature":
            hander = FeatureHandler()
        elif label == "psm":
            hander = PSMHandler()
        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
            # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()


def read_large_parquet(parquet_path: str, batch_size: int = 500000):
    """_summary_
    :param parquet_path: _description_
    :param batch_size: _description_, defaults to 100000
    :yield: _description_
    """
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batch_df = batch.to_pandas()
        yield batch_df


# register_file
def register_file_to_json(project_file, attach_file, category, replace_existing):
    register = ProjectHandler(project_json_file=project_file)
    register.add_quantms_file(attach_file, category, replace_existing)
    register.save_updated_project_info(output_file_name=project_file)


# get best_scan_number
def load_best_scan_number(diann_psm_path: str, diann_feature_path: str, output_path: str):
    p = Parquet(diann_psm_path)
    psm_df = p.load_psm_scan()
    pqwriter = None
    hander = FeatureHandler()
    for df in read_large_parquet(diann_feature_path):
        df = df.merge(psm_df, on=["peptidoform", "charge"], how="left")
        df.drop(["best_psm_scan_number"], inplace=True, axis=1)
        df.rename(columns={"scan_number": "best_psm_scan_number"}, inplace=True)

        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
            # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()


def fill_start_and_end(row, protein_dict):
    start = []
    end = []
    for protein in list(row["protein_accessions"]):
        automaton = ahocorasick.Automaton()
        automaton.add_word(row["sequence"], row["sequence"])
        automaton.make_automaton()
        for item in automaton.iter(protein_dict[protein]):
            end.append(item[0])
            start.append(item[0] - len(row["sequence"]) + 1)
    return start, end


def generate_start_and_end_from_fasta(parquet_path, fasta_path, label, output_path):
    p = Parquet(parquet_path)
    protein_dict = p.get_protein_dict(fasta_path)
    if label == "feature":
        hander = FeatureHandler()
    elif label == "psm":
        hander = PSMHandler()
    pqwriter = None
    for df in read_large_parquet(parquet_path):
        df[["protein_start_positions", "protein_end_positions"]] = df[["sequence", "protein_accessions"]].swifter.apply(
            lambda row: fill_start_and_end(row, protein_dict),
            axis=1,
            result_type="expand",
        )
        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
            # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()


# report


def convert_to_base64(fig):
    figfile = BytesIO()
    fig.figure.savefig(figfile, format="png")
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue())
    figdata_str = str(figdata_png, "utf-8")

    return "data:image/png;base64," + figdata_str


def get_msgs_from_project(path):
    f = open(path, "r")
    content = f.read()
    res = json.loads(content)
    return (
        res["project_accession"],
        res["project_title"],
        res["project_data_description"],
    )


def generate_project_report(project_folder):
    msgs = {
        "project": "",
        "projectTitle": "",
        "projectDescription": "",
        "featureProtrins": "",
        "featurePeptides": "",
        "featureSamples": "",
        "featurePeptidoforms": "",
        "featureMsruns": "",
        "featureImg1": "",
        "featureImg2": "",
        "featureImg3": "",
        "psmProtrins": "",
        "psmPeptides": "",
        "psmPeptidoforms": "",
        "psms": "",
        "psmMsruns": "",
        "psmImg1": "",
        "aeProtrins": "",
        "aeSamples": "",
        "aeImg1": "",
    }
    file_list = os.listdir(project_folder)
    os.chdir(project_folder)
    if "project.json" in file_list:
        project_info = get_msgs_from_project("project.json")
        msgs["project"] = project_info[0]
        msgs["projectTitle"] = project_info[1]
        msgs["projectDescription"] = project_info[2]
    feature_paths = [f for f in file_list if f.endswith(".feature.parquet")]
    if len(feature_paths) > 0:
        feature_statistics = ParquetStatistics(feature_paths[0])
        msgs["featureProtrins"] = feature_statistics.get_number_of_proteins()
        msgs["featurePeptides"] = feature_statistics.get_number_of_peptides()
        msgs["featureSamples"] = feature_statistics.get_number_of_samples()
        msgs["featurePeptidoforms"] = feature_statistics.get_number_of_peptidoforms()
        msgs["featureMsruns"] = feature_statistics.get_number_msruns()
        msgs["featureImg1"] = convert_to_base64(plot_intensity_distribution_of_samples(feature_paths[0]))
        msgs["featureImg2"] = convert_to_base64(plot_peptide_distribution_of_protein(feature_paths[0]))
        msgs["featureImg3"] = convert_to_base64(plot_intensity_box_of_samples(feature_paths[0]))
    psm_paths = [f for f in file_list if f.endswith(".psm.parquet")]
    if len(psm_paths) > 0:
        psm_statistics = ParquetStatistics(psm_paths[0])
        msgs["psmProtrins"] = psm_statistics.get_number_of_proteins()
        msgs["psmPeptides"] = psm_statistics.get_number_of_peptides()
        msgs["psmPeptidoforms"] = psm_statistics.get_number_of_peptidoforms()
        msgs["psms"] = psm_statistics.get_number_of_psms()
        msgs["psmMsruns"] = psm_statistics.get_number_msruns()
        sdrf_paths = [f for f in file_list if f.endswith(".sdrf.tsv")]
        if len(sdrf_paths) > 0:
            msgs["psmImg1"] = convert_to_base64(plot_peptides_of_lfq_condition(psm_paths[0], sdrf_paths[0]))
    ae_paths = [f for f in file_list if f.endswith(".absolute.tsv")]
    if len(ae_paths) > 0:
        absolute_stats = IbaqStatistics(ibaq_path=ae_paths[0])
        msgs["aeProtrins"] = absolute_stats.get_number_of_proteins()
        msgs["aeSamples"] = absolute_stats.get_number_of_samples()
        msgs["aeImg1"] = convert_to_base64(plot_distribution_of_ibaq(ae_paths[0]))

    template_str = string.Template(report)
    with open("report.html", "w", encoding="utf8") as f:
        f.write(template_str.substitute(msgs))
