import os
import re

from collections import defaultdict

import ahocorasick
import duckdb

import mygene
import pandas as pd
import pyarrow.parquet as pq
from Bio import SeqIO

from quantmsio.core.openms import OpenMSHandler

from quantmsio.utils.pride_utils import generate_gene_name_map
from quantmsio.utils.pride_utils import get_gene_accessions
from quantmsio.utils.pride_utils import get_unanimous_name


def check_string(re_exp, strings):
    res = re.search(re_exp, strings)
    if res:
        return True
    else:
        return False


def map_spectrum_mz(mz_path: str, scan: str, mzml: dict, mzml_directory: str):
    """
    mz_path: mzML file path
    scan: scan number
    mzml: OpenMSHandler object
    """
    reference = mz_path
    mz_path = os.path.join(mzml_directory, mz_path + ".mzML")
    number_peaks, mz_array, intensity_array = mzml[reference].get_spectrum_from_scan(mz_path, int(scan))
    return number_peaks, mz_array, intensity_array


def fill_start_and_end(row, protein_dict):
    """
    Map seq location from fasta file.
    return: Tuple of location
    """
    positions = []
    for protein in list(row["pg_accessions"]):
        automaton = ahocorasick.Automaton()
        automaton.add_word(row["sequence"], row["sequence"])
        automaton.make_automaton()
        if protein not in protein_dict:
            return None
        for item in automaton.iter(protein_dict[protein]):
            end = item[0]
            start = item[0] - len(row["sequence"]) + 1
            positions.append(f"{start}:{end}")
    return positions


class Query:

    def __init__(self, parquet_path: str):
        if os.path.exists(parquet_path):
            self._path = parquet_path
            self.parquet_db = duckdb.connect(config={"max_memory": "16GB", "worker_threads": 4})
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path)
            )
        else:
            raise FileNotFoundError(f"the file {parquet_path} does not exist.")

    def get_report_from_database(self, runs: list, columns: list = None):
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        cols = cols.replace("unique", '"unique"')
        database = self.parquet_db.sql(
            """
            select {} from parquet_db
            where reference_file_name IN {}
            """.format(
                cols, tuple(runs)
            )
        )
        report = database.df()
        return report

    def get_samples_from_database(self, samples: list, columns: list = None):
        """
        This function loads the report from the duckdb database for a group of samples.
        :param runs: A list of samples
        :return: The report
        """
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        cols = cols.replace("unique", '"unique"')
        database = self.parquet_db.sql(
            """
            select {} from parquet_db
            where sample_accession IN {}
            """.format(
                cols, tuple(samples)
            )
        )
        report = database.df()
        return report

    def iter_samples(self, file_num: int = 20, columns: list = None):
        """
        :params file_num: The number of files being processed at the same time(default 10)
        :yield: _description_
        """
        samples = self.get_unique_samples()
        ref_list = [samples[i : i + file_num] for i in range(0, len(samples), file_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs, columns)
            yield refs, batch_df

    def iter_chunk(self, batch_size: int = 500000, columns: list = None):
        """_summary_
        :param batch_size: _description_, defaults to 100000
        :yield: _description_
        """
        parquet_file = pq.ParquetFile(self._path)
        for batch in parquet_file.iter_batches(batch_size=batch_size, columns=columns):
            batch_df = batch.to_pandas()
            yield batch_df

    def iter_file(self, file_num: int = 10, columns: list = None):
        """
        :params file_num: The number of files being processed at the same time(default 10)
        :yield: _description_
        """
        references = self.get_unique_references()
        ref_list = [references[i : i + file_num] for i in range(0, len(references), file_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs, columns)
            yield refs, batch_df

    def get_spectrum_msg(self, reference: str, scan: str, mzml_directory: str):
        """
        :params reference: reference_file_name
        :params scan: scan
        :params mzml_directory: Mzml folder
        :return (number_peaks, mz_array, intensity_array)
        """
        mzml_handler = {reference: OpenMSHandler()}
        return map_spectrum_mz(reference, scan, mzml_handler, mzml_directory)

    def inject_position_msg(self, df: pd.DataFrame, protein_dict: dict):
        """
        :params df: parquet file
        :params protein_dict: {protein_accession:seq}
        :retrun df
        """
        df["pg_positions"] = df[["sequence", "pg_accessions"]].apply(
            lambda row: fill_start_and_end(row, protein_dict), axis=1
        )
        return df

    def inject_gene_msg(
        self,
        df: pd.DataFrame,
        map_gene_names: dict,
        species: str = "human",
    ):
        """
        :params df: parquet file
        :params fasta: refence fasta file
        :params map_parameter: map_protein_name or map_protein_accession
        :params species: default human
        :return df
        """
        if "pg_accessions" in df.columns:
            df["gg_names"] = df["pg_accessions"].apply(lambda x: get_unanimous_name(x, map_gene_names))
        else:
            df["gg_names"] = df["mp_accessions"].apply(lambda x: get_unanimous_name(x, map_gene_names))
        gene_list = list(set([gene for gene_names in df["gg_names"] if gene_names is not None for gene in gene_names]))
        gene_accessions = self.get_gene_accessions(gene_list, species)
        df["gg_accessions"] = df["gg_names"].apply(lambda x: get_gene_accessions(x, gene_accessions))

        return df

    def get_protein_to_gene_map(self, fasta: str, map_parameter: str = "map_protein_accession"):
        map_gene_names = generate_gene_name_map(fasta, map_parameter)
        return map_gene_names

    def get_protein_dict(self, fasta_path):
        """
        return: protein_map {protein_accession:seq}
        """
        df = self.parquet_db.sql("SELECT DISTINCT pg_accessions FROM parquet_db").df()
        proteins = set()
        for protein_accessions in df["pg_accessions"].tolist():
            proteins.update(set(protein_accessions))
        protein_dict = {}
        for seq in SeqIO.parse(fasta_path, "fasta"):
            p_name = seq.id.split("|")[1]
            if p_name in proteins and p_name not in protein_dict:
                protein_dict[p_name] = str(seq.seq)
        return protein_dict

    def load_psm_scan(self):
        psm_df = self.parquet_db.sql(
            """
            SELECT peptidoform,charge,scan_number
            FROM (
            SELECT peptidoform,charge,scan_number, ROW_NUMBER()
            OVER (PARTITION BY peptidoform,charge ORDER BY global_qvalue ASC) AS row_num
            FROM parquet_db
            ) AS subquery
            WHERE row_num = 1;
            """
        ).df()
        return psm_df

    def get_unique_references(self):
        """
        return: A list of deduplicated reference
        """
        unique_reference = self.parquet_db.sql("SELECT DISTINCT reference_file_name FROM parquet_db").df()

        return unique_reference["reference_file_name"].tolist()

    def get_unique_peptides(self):
        """
        return: A list of deduplicated peptides.
        """
        unique_peps = self.parquet_db.sql("SELECT DISTINCT sequence FROM parquet_db").df()

        return unique_peps["sequence"].tolist()

    def get_unique_proteins(self):
        """
        return: A list of deduplicated proteins.
        """

        unique_prts = self.parquet_db.sql("SELECT mp_accessions FROM parquet_db").df()
        proteins = set()
        for protein_accessions in unique_prts["mp_accessions"]:
            if protein_accessions is not None:
                proteins.update(set(protein_accessions))
        return list(proteins)

    def get_unique_genes(self):
        """
        return: A list of deduplicated genes.
        """

        unique_prts = self.parquet_db.sql("SELECT DISTINCT gg_names FROM parquet_db").df()

        return unique_prts["gg_names"].tolist()

    def get_unique_samples(self):
        """
        return: A list of deduplicated sampless.
        """
        unique_peps = self.parquet_db.sql("SELECT DISTINCT sample_accession FROM parquet_db").df()
        return unique_peps["sample_accession"].tolist()

    def query_peptide(self, peptide: str, columns: list = None):
        """
        peptide: Peptide that need to be queried.
        return: A DataFrame of all information about query peptide.
        """

        if check_string("^[A-Z]+$", peptide):
            cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
            return self.parquet_db.sql(f"SELECT {cols} FROM parquet_db WHERE sequence ='{peptide}'").df()
        else:
            raise KeyError("Illegal peptide!")

    def query_peptides(self, peptides: list, columns: list = None):
        """
        :params protein: Protein that need to be queried.
        return: A DataFrame of all information about query proteins.
        """
        for p in peptides:
            if not check_string("^[A-Z]+$", p):
                raise KeyError("Illegal peptide!")
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        database = self.parquet_db.sql(f"select {cols} from parquet_db where sequence IN {tuple(peptides)}")
        return database.df()

    def query_proteins(self, proteins: list, columns: list = None):
        """
        :params protein: Protein that need to be queried.
        return: A DataFrame of all information about query proteins.
        """
        for p in proteins:
            if not check_string("^[A-Z]+", p):
                raise KeyError("Illegal protein!")
        proteins_key = [f"pg_accessions LIKE '%{p}%'" for p in proteins]
        query_key = " OR ".join(proteins_key)
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        database = self.parquet_db.sql(f"SELECT {cols} FROM parquet_db WHERE {query_key}")
        return database.df()

    def query_protein(self, protein: str, columns: list = None):
        """
        :params protein: Protein that need to be queried.
        return: A DataFrame of all information about query protein.
        """
        cols = ", ".join(columns) if columns and isinstance(columns, list) else "*"
        if check_string("^[A-Z]+", protein):
            return self.parquet_db.sql(f"SELECT {cols} FROM parquet_db WHERE pg_accessions LIKE '%{protein}%'").df()
        else:
            raise KeyError("Illegal protein!")

    def get_gene_accessions(self, gene_list: list, species: str = "human"):
        """
        :params gene_list
        """
        mg = mygene.MyGeneInfo()
        gene_accessions = mg.querymany(gene_list, scopes="symbol", species=species, fields="accession")
        gene_accessions_maps = defaultdict(list)
        for obj in gene_accessions:
            if "accession" in obj and "genomic" in obj["accession"]:
                gene_accessions_maps[obj["query"]] = ",".join(obj["accession"]["genomic"])
        return gene_accessions_maps
