"""
tools.py provide some function for generate optional feature:
mz_array: 
intensity_array:
num_peaks:
gene_names:
gene_accessions:
compare parquet
plot veen
"""
import os
import re

import pandas as pd
import matplotlib.pyplot as plt
from venn import venn

import pyarrow as pa
import pyarrow.parquet as pq
from quantms_io.core.feature import FeatureHandler
from quantms_io.core.openms import OpenMSHandler



def extract_len(fle, header):
    map_tag = {"PSH": "PSM", "PEH": "PEP", "PRH": "PRT"}
    if os.stat(fle).st_size == 0:
        raise ValueError("File is empty")
    f = open(fle)
    pos = 0
    line = f.readline()
    while line.split("\t")[0] != header:
        pos = f.tell()
        line = f.readline()
    line = f.readline()
    fle_len = 0
    while line.split("\t")[0] == map_tag[header]:
        fle_len += 1
        line = f.readline()
    f.close()
    return fle_len, pos


# optional about spectrum
def map_spectrum_mz(mz_path: str, scan: str, Mzml: OpenMSHandler, mzml_directory: str):
    """
    mz_path: mzML file path
    scan: scan number
    mzml: OpenMSHandler object
    """
    if mzml_directory.endswith("/"):
        mz_path = mzml_directory + mz_path + ".mzML"
    else:
        mz_path = mzml_directory + "/" + mz_path + ".mzML"
    mz_array, array_intensity = Mzml.get_spectrum_from_scan(mz_path, int(scan))
    return mz_array, array_intensity, 0


def generate_features_of_spectrum(parquet_path: str, mzml_directory: str,output_path:str):
    '''
    parquet_path: parquet file path
    mzml_directory: mzml file directory path
    '''
    table = pd.read_parquet(parquet_path)
    Mzml = OpenMSHandler()

    table[["mz", "array_intensity", "num_peaks"]] = table[
        ["best_psm_reference_file_name", "best_psm_scan_number"]
    ].apply(
        lambda x: map_spectrum_mz(
            x["best_psm_reference_file_name"],
            x["best_psm_scan_number"],
            Mzml,
            mzml_directory,
        ),
        axis=1,
        result_type="expand",
    )
    Feature = FeatureHandler()
    parquet_table = pa.Table.from_pandas(table, schema=Feature.schema)
    pq.write_table(parquet_table, output_path, compression="snappy", store_schema=True)


# option about gene
def generate_gene_names(parquet_path: str, mzTab_path: str):
    """
    get gene names from mzTab into parquet.
    """
    gene_map = extract_optional_gene_dict(mzTab_path)
    table = pd.read_parquet(parquet_path)
    keys = gene_map.keys()
    table["gene_names"] = table["protein_accessions"].apply(
        lambda x: gene_map[",".join(x)] if ",".join(x) in keys else pd.NA
    )
    Feature = FeatureHandler()
    parquet_table = pa.Table.from_pandas(table, schema=Feature.schema)

    pq.write_table(parquet_table, parquet_path, compression="snappy", store_schema=True)


def query_accessions(gene_name, Email, species: str = "Homo sapiens", OX: str = "9606"):
    """
    Getting the gene_accessions is time-consuming work, so we support single gene query.
    gene_name: 'gene name'
    Email: 'a available email address'
    species: ''
    OX: 'Organism Taxonomy ID'
    """
    from Bio import Entrez
    Entrez.email = Email
    q_name = gene_name + " AND " + species + "[porgn:__txid" + OX + "]"
    handle = Entrez.esearch(db="nucleotide", term=q_name)
    record = Entrez.read(handle)
    id_list = record["IdList"]
    accessions = get_accessions(id_list)
    return accessions


def get_accessions(id_list):
    from Bio import Entrez, SeqIO
    accessions = []
    for id_q in id_list:
        handle = Entrez.efetch(db="nucleotide", id=id_q, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        accessions.append(record.id)
    return accessions


def load_prt(mzTab_path, **kwargs):
    fle_len, pos = extract_len(mzTab_path, "PRH")
    f = open(mzTab_path)
    f.seek(pos)
    return pd.read_csv(f, nrows=fle_len, **kwargs)


def extract_optional_gene_dict(mzTab_path):
    """
    return: a dict, key is protein accessions, value is gene names
    """
    PRT = load_prt(mzTab_path, "PRH", sep="\t", usecols=["accession", "description"])
    PRT.loc[:, "gene_names"] = (
        PRT["description"]
        .fillna("unknown")
        .apply(
            lambda x: re.findall(r"GN=(\w+)", x)
            if re.findall(r"GN=(\w+)", x)
            else pd.NA
        )
    )
    gene_map = (
        PRT[["accession", "gene_names"]]
        .set_index("accession")
        .to_dict(orient="dict")["gene_names"]
    )

    return gene_map

#plot venn
def plot_peptidoform_charge_venn(parquet_path_list,labels):
    data_map = {}
    for parquet_path,label in zip(parquet_path_list,labels):
        df = pd.read_parquet(parquet_path,columns=['peptidoform','charge'])
        psm_message = 'Total number of PSM for ' + label + ': ' + str(len(df))
        print(psm_message)
        unique_pep_forms = set((df['peptidoform'] + df['charge'].astype(str)).to_list())
        pep_form_message = 'Total number of Peptidoform for ' + label + ': ' + str(len(unique_pep_forms))
        print(pep_form_message)
        data_map[label] = unique_pep_forms
    plt.figure(figsize=(16, 12), dpi=500)
    venn(data_map, legend_loc="upper right",figsize=(16, 12))
    plt.savefig('pep_form_compare_venn.png')
    plt.show()

def plot_sequence_venn(parquet_path_list,labels):
    data_map = {}
    for parquet_path,label in zip(parquet_path_list,labels):
        sequence = pd.read_parquet(parquet_path,columns=['sequence'])
        unique_seqs = set(sequence['sequence'].to_list())
        pep_message = 'Total number of peptide for ' + label + ": " + str(unique_seqs)
        print(pep_message)
        data_map[label] = unique_seqs
    plt.figure(figsize=(16, 12), dpi=500)
    venn(data_map, legend_loc="upper right",figsize=(16, 12))
    plt.savefig('sequence_compare_venn.png')
    plt.show()
    