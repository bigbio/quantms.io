"""
tools.py provide some function for generate optional feature:
mz_array: 
intensity_array:
num_peaks:
gene_names:
gene_accessions:
compare parquet
plot veen
...
"""
import os
import re
import json

import pandas as pd
import matplotlib.pyplot as plt
from venn import venn

import pyarrow as pa
import pyarrow.parquet as pq
from quantms_io.core.feature import FeatureHandler
from quantms_io.core.openms import OpenMSHandler
from quantms_io.core.psm import PSMHandler
from quantms_io.core.project import ProjectHandler
from quantms_io.utils.file_utils import extract_len
from quantms_io.core.project import create_uuid_filename
import duckdb


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


def generate_features_of_spectrum(parquet_path: str, mzml_directory: str,output_path:str,label,chunksize,partition:str=None):
    """
    parquet_path: parquet file path
    mzml_directory: mzml file directory path
    """
    pqwriters = {}
    pqwriter_no_part = None
    for table in read_large_parquet(parquet_path,batch_size=chunksize):
        #table = pd.read_parquet(parquet_path)
        mzml = OpenMSHandler()
        if label == 'feature':
            table[["mz_array", "intensity_array", "num_peaks"]] = table[
                ["best_psm_reference_file_name", "best_psm_scan_number"]
            ].apply(
                lambda x: map_spectrum_mz(
                    x["best_psm_reference_file_name"],
                    x["best_psm_scan_number"],
                    mzml,
                    mzml_directory,
                ),
                axis=1,
                result_type="expand",
            )
            hander = FeatureHandler()
        else:
            table[["mz_array", "intensity_array", "num_peaks"]] = table[
                ["reference_file_name", "scan_number"]
            ].apply(
                lambda x: map_spectrum_mz(
                    x["reference_file_name"],
                    x["scan_number"],
                    mzml,
                    mzml_directory,
                ),
                axis=1,
                result_type="expand",
            )
            hander = PSMHandler()
        #save
        if partition == 'charge':
            for key, df in table.groupby(['charge']):
                parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
                save_path = output_path.split(".")[0] +"-" + str(key) + '.parquet'
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        elif partition == 'reference_file_name':
            for key, df in table.groupby(['reference_file_name']):
                parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
                save_path = output_path.split(".")[0] + "-" + str(key) + '.parquet'
                if not os.path.exists(save_path):
                    pqwriter = pq.ParquetWriter(save_path, parquet_table.schema)
                    pqwriters[key] = pqwriter
                pqwriters[key].write_table(parquet_table)
        else:
            parquet_table = pa.Table.from_pandas(table, schema=hander.schema)
            if not pqwriter_no_part:
                pqwriter_no_part = pq.ParquetWriter(output_path, parquet_table.schema)
            pqwriter_no_part.write_table(parquet_table)
    #close f
    if not partition:
        if pqwriter_no_part:
            pqwriter_no_part.close()
    else:
        for pqwriter in pqwriters.values():
            pqwriter.close()


# option about gene
def generate_gene_names(parquet_path: str, mz_tab_path: str):
    """
    get gene names from mzTab into parquet.
    """
    gene_map = extract_optional_gene_dict(mz_tab_path)
    table = pd.read_parquet(parquet_path)
    keys = gene_map.keys()
    table["gene_names"] = table["protein_accessions"].apply(
        lambda x: gene_map[",".join(x)] if ",".join(x) in keys else pd.NA
    )
    feature = FeatureHandler()
    parquet_table = pa.Table.from_pandas(table, schema=feature.schema)

    pq.write_table(parquet_table, parquet_path, compression="snappy", store_schema=True)


def query_accessions(gene_name, email, species: str = "Homo sapiens", ox: str = "9606"):
    """
    getting the gene_accessions is time-consuming work, so we support single gene query.
    gene_name: 'gene name'
    Email: 'a available email address'
    species: ''
    OX: 'Organism Taxonomy ID'
    """
    from Bio import Entrez
    Entrez.email = email
    q_name = gene_name + " AND " + species + "[porgn:__txid" + ox + "]"
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


def load_prt(mz_tab_path, **kwargs):
    fle_len, pos = extract_len(mz_tab_path, "PRH")
    f = open(mz_tab_path)
    f.seek(pos)
    return pd.read_csv(f, nrows=fle_len, **kwargs)


def extract_optional_gene_dict(mz_tab_path):
    """
    return: a dict, key is protein accessions, value is gene names
    """
    prt = load_prt(mz_tab_path, "PRH", sep="\t", usecols=["accession", "description"])
    prt.loc[:, "gene_names"] = (
        prt["description"]
        .fillna("unknown")
        .apply(
            lambda x: re.findall(r"GN=(\w+)", x)
            if re.findall(r"GN=(\w+)", x)
            else pd.NA
        )
    )
    gene_map = (
        prt[["accession", "gene_names"]]
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
    venn(data_map, legend_loc="upper right",figsize=(16, 12),fmt="{size}({percentage:.1f}%)")
    plt.savefig('pep_form_compare_venn.png')

def plot_sequence_venn(parquet_path_list,labels):
    data_map = {}
    for parquet_path,label in zip(parquet_path_list,labels):
        sequence = pd.read_parquet(parquet_path,columns=['sequence'])
        unique_seqs = set(sequence['sequence'].to_list())
        pep_message = 'Total number of peptide for ' + label + ": " + str(len(unique_seqs))
        print(pep_message)
        data_map[label] = unique_seqs
    plt.figure(figsize=(16, 12), dpi=500)
    venn(data_map, legend_loc="upper right",figsize=(16, 12),fmt="{size}({percentage:.1f}%)")
    plt.savefig('sequence_compare_venn.png')

# gei unqnimous name
from collections import defaultdict
def map_protein_for_parquet(parquet_path,fasta,output_path,map_parameter,label):
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
    if map_parameter == 'map_protein_name':
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession,name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(name)
            map_protein_names[accession].add(name)
            map_protein_names[name].add(name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession,name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(accession)
            map_protein_names[accession].add(accession)
            map_protein_names[name].add(accession)
    change_and_save_parquet(parquet_path,map_protein_names,output_path,label)

def map_protein_for_tsv(path: str,fasta: str, output_path: str, map_parameter: str):
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param path: de_path or ae_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    from Bio import SeqIO
    map_protein_names = defaultdict(set)
    if map_parameter == 'map_protein_name':
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession,name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(name)
            map_protein_names[accession].add(name)
            map_protein_names[name].add(name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession,name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(accession)
            map_protein_names[accession].add(accession)
            map_protein_names[name].add(accession)
    df,content = load_de_or_ae(path)
    df['protein'] = df['protein'].apply(lambda x: get_unanimous_name(x,map_protein_names))
    content += df.columns.str.cat(sep="\t") + "\n"
    for index, row in df.iterrows():
        content += '\t'.join(map(str, row)).strip() + "\n"
    with open(output_path, "w",encoding='utf8') as f:
        f.write(content)

def load_de_or_ae(path):
    f = open(path,encoding='utf-8')
    line = f.readline()
    pos = 0 
    content = ""
    while line.startswith('#'):
        pos = f.tell()
        content += line
        line = f.readline()
    f.seek(pos-1)
    return pd.read_csv(f,sep='\t'),content
        
def change_and_save_parquet(parquet_path,map_dict,output_path,label):
    pqwriter = None
    for df in read_large_parquet(parquet_path):
        df['protein_accessions'] = df['protein_accessions'].apply(lambda x: get_unanimous_name(x,map_dict))
        
        if label == 'feature':
            hander = FeatureHandler()
        elif label == 'psm':
            hander = PSMHandler()
        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
        # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close()
        
def get_unanimous_name(protein_accessions,map_dict):
    if isinstance(protein_accessions,str):
        if ';' in protein_accessions:
            protein_accessions = protein_accessions.split(";")
        else:
            protein_accessions = protein_accessions.split(",")
    unqnimous_names = []
    for accession in protein_accessions:
        unqnimous_names.append(list(map_dict[accession])[0])
    return unqnimous_names
        
def read_large_parquet(parquet_path: str, batch_size: int = 500000):
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batch_df = batch.to_pandas()
        yield batch_df

#register_file
def register_file_to_json(project_file,attach_file,category,replace_existing):
    register= ProjectHandler(project_json_file=project_file)
    register.add_quantms_file(attach_file,category,replace_existing)
    register.save_updated_project_info(output_file_name=project_file)

#check result of psms or features
def generate_report_of_psms_or_features(check_dir,label):
    if not os.path.exists(check_dir):
        raise Exception("not file path")
    file_list = os.listdir(check_dir)
    if label == 'psm':
        check_list = [file for file in file_list if file.endswith(".psm.parquet")]
    elif label == 'feature':
        check_list = [file for file in file_list if file.endswith(".feature.parquet")]
    
    output_lines = ''
    for file_path in check_list:
        output_lines += 'Name: ' + file_path + '\n'
        file_path = check_dir + "/" + file_path
        file_size = get_file_size(file_path)
        output_lines += 'File size: ' + file_size + '\n'
        df = pd.read_parquet(file_path,columns=['protein_accessions','peptidoform','charge'])
        output_lines += 'Total number of Peptides: ' + str(len(df.groupby(['peptidoform','charge']))) + '\n'
        proteins = set()
        df['protein_accessions'].apply(lambda x: proteins.update(set(x)))
        output_lines += 'Total number of Proteins: ' + str(len(proteins)) + '\n\n'

    output_path =  create_uuid_filename(label+'s_report','.txt')
    with open(output_path, "w",encoding='utf8') as f:
            f.write(output_lines)

def get_file_size(file_path):
    fsize = os.path.getsize(file_path)
    if fsize < 1024:
        return str(round(fsize,2)) + 'Byte'
    else: 
        kbx = fsize/1024
        if kbx < 1024:
            return str(round(kbx,2)) + 'K'
        else:
            mbx = kbx /1024
            if mbx < 1024:
                return str(round(mbx,2)) + 'M'
            else:
                return str(round(mbx/1024)) + 'G'
                
#covert ae or de to json
def convert_to_json(file_path):
    """
    by providing the json format of AE and DE files for retrieval. return json
    """
    table,content = load_de_or_ae(file_path)
    output = {}
    pattern = r'[\\|\|//|/]'
    file_name = re.split(pattern,file_path)[-1]
    output['id'] = file_name
    output['metadata'] = content
    records = {}
    for col in table.columns:
        records[col] = table.loc[:,col].to_list()
    output['records'] = records
    b = json.dumps(output)
    output_path = ".".join(file_name.split('.')[:-1]) + '.json'
    f = open(output_path, 'w')
    f.write(b)
    f.close()

#get best_scan_number
def load_best_scan_number(diann_psm_path:str,diann_feature_path:str,output_path:str):
    database = duckdb.connect(config={
    "max_memory" : "16GB",
    "worker_threads": 4
    })
    psm_table = database.read_parquet(diann_psm_path)
    psm_df = database.execute(
        """    
        SELECT peptidoform,charge,scan_number
        FROM (
        SELECT peptidoform,charge,scan_number, ROW_NUMBER() OVER (PARTITION BY peptidoform,charge ORDER BY global_qvalue ASC) AS row_num
        FROM psm_table
        ) AS subquery
        WHERE row_num = 1;
        """
    ).df()
    pqwriter = None
    hander = FeatureHandler()
    for df in read_large_parquet(diann_feature_path):
        df = df.merge(psm_df,on=['peptidoform','charge'],how='left')
        df.drop(["best_psm_scan_number"],inplace=True,axis=1)
        df.rename(columns={
            "scan_number": "best_psm_scan_number"
        },inplace=True)

        parquet_table = pa.Table.from_pandas(df, schema=hander.schema)
        if not pqwriter:
        # create a parquet write object giving it an output file
            pqwriter = pq.ParquetWriter(output_path, parquet_table.schema)
        pqwriter.write_table(parquet_table)
    if pqwriter:
        pqwriter.close() 

