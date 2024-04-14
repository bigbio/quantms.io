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
    if mzml_directory.endswith("/"):
        mz_path = mzml_directory + mz_path + ".mzML"
    else:
        mz_path = mzml_directory + "/" + mz_path + ".mzML"
    mz_array, array_intensity = mzml[reference].get_spectrum_from_scan(mz_path, int(scan))
    return mz_array, array_intensity, 0


def fill_start_and_end(row, protein_dict):
    """
    Map seq location from fasta file.
    return: Tuple of location
    """
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


class Parquet:

    def __init__(self, parquet_path: str):
        if os.path.exists(parquet_path):
            self._path = parquet_path
            self.parquet_db = duckdb.connect(config={"max_memory": "16GB", "worker_threads": 4})
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path)
            )
        else:
            raise FileNotFoundError(f"the file {parquet_path} does not exist.")

    def get_report_from_database(self, runs: list):
        """
        This function loads the report from the duckdb database for a group of ms_runs.
        :param runs: A list of ms_runs
        :return: The report
        """
        database = self.parquet_db.sql(
            """
            select * from parquet_db
            where reference_file_name IN {}
            """.format(
                tuple(runs)
            )
        )
        report = database.df()
        return report

    def get_samples_from_database(self, samples: list):
        """
        This function loads the report from the duckdb database for a group of samples.
        :param runs: A list of samples
        :return: The report
        """
        database = self.parquet_db.sql(
            """
            select * from parquet_db
            where sample_accession IN {}
            """.format(
                tuple(samples)
            )
        )
        report = database.df()
        return report

    def iter_samples(self, file_num: int = 20):
        """
        :params file_num: The number of files being processed at the same time(default 10)
        :yield: _description_
        """
        samples = self.get_unique_samples()
        ref_list = [samples[i : i + file_num] for i in range(0, len(samples), file_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs)
            yield refs, batch_df

    def iter_chunk(self, batch_size: int = 500000):
        """_summary_
        :param batch_size: _description_, defaults to 100000
        :yield: _description_
        """
        parquet_file = pq.ParquetFile(self._path)
        for batch in parquet_file.iter_batches(batch_size=batch_size):
            batch_df = batch.to_pandas()
            yield batch_df

    def iter_file(self, file_num: int = 10):
        """
        :params file_num: The number of files being processed at the same time(default 10)
        :yield: _description_
        """
        references = self.get_unique_references()
        ref_list = [references[i : i + file_num] for i in range(0, len(references), file_num)]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs)
            yield batch_df

    def inject_spectrum_msg(self, df: pd.DataFrame, mzml_directory: str):
        """
        :params df: parquet file
        :params maml_directory: Mzml folder
        :return df
        """
        refs = df["reference_file_name"].unique()
        mzml = {ref: OpenMSHandler() for ref in refs}
        if "best_psm_reference_file_name" in df.columns:
            df[["mz_array", "intensity_array", "num_peaks"]] = df[
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
            return df
        elif "reference_file_name" in df.columns:
            df[["mz_array", "intensity_array", "num_peaks"]] = df[["reference_file_name", "scan_number"]].apply(
                lambda x: map_spectrum_mz(
                    x["reference_file_name"],
                    x["scan_number"],
                    mzml,
                    mzml_directory,
                ),
                axis=1,
                result_type="expand",
            )
            return df

    def inject_position_msg(self, df: pd.DataFrame, protein_dict: dict):
        """
        :params df: parquet file
        :params protein_dict: {protein_accession:seq}
        :retrun df
        """
        df[["protein_start_positions", "protein_end_positions"]] = df[["sequence", "protein_accessions"]].apply(
            lambda row: fill_start_and_end(row, protein_dict),
            axis=1,
            result_type="expand",
        )
        return df

    def inject_gene_msg(
        self,
        df: pd.DataFrame,
        fasta: str,
        map_parameter: str = "map_protein_accession",
        species: str = "human",
    ):
        """
        :params df: parquet file
        :params fasta: refence fasta file
        :params map_parameter: map_protein_name or map_protein_accession
        :params species: default human
        :return df
        """
        map_gene_names = generate_gene_name_map(fasta, map_parameter)
        df["gene_names"] = df["protein_accessions"].apply(lambda x: get_unanimous_name(x, map_gene_names))
        gene_list = self.get_gene_list(map_gene_names)
        gene_accessions = self.get_gene_accessions(gene_list, species)
        df["gene_accessions"] = df["gene_names"].apply(lambda x: get_gene_accessions(x, gene_accessions))

        return df

    def get_protein_dict(self, fasta_path):
        """
        return: protein_map {protein_accession:seq}
        """
        df = self.parquet_db.sql("SELECT DISTINCT protein_accessions FROM parquet_db").df()
        proteins = set()
        for protein_accessions in df["protein_accessions"].tolist():
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

        unique_prts = self.parquet_db.sql("SELECT DISTINCT protein_accessions FROM parquet_db").df()

        return unique_prts["protein_accessions"].tolist()

    def get_unique_genes(self):
        """
        return: A list of deduplicated genes.
        """

        unique_prts = self.parquet_db.sql("SELECT DISTINCT gene_names FROM parquet_db").df()

        return unique_prts["gene_names"].tolist()

    def get_unique_samples(self):
        """
        return: A list of deduplicated sampless.
        """
        unique_peps = self.parquet_db.sql("SELECT DISTINCT sample_accession FROM parquet_db").df()
        return unique_peps["sample_accession"].tolist()

    def query_peptide(self, peptide: str):
        """
        peptide: Peptide that need to be queried.
        return: A DataFrame of all information about query peptide.
        """

        if check_string("^[A-Z]+$", peptide):
            return self.parquet_db.sql(f"SELECT * FROM parquet_db WHERE sequence ='{peptide}'").df()
        else:
            return KeyError("Illegal peptide!")

    def query_protein(self, protein: str):
        """
        :params protein: Protein that need to be queried.
        return: A DataFrame of all information about query protein.
        """
        if check_string("^[A-Z]+", protein):
            return self.parquet_db.sql(f"SELECT * FROM parquet_db WHERE protein_accessions ='{protein}'").df()
        else:
            return KeyError("Illegal protein!")

    def get_gene_list(self, map_gene_names: dict):
        """
        :params map_gene_names: protenin => gene
        return: unique gene list
        """
        unique_prts = self.get_unique_proteins()
        gene_names = [get_unanimous_name(proteins, map_gene_names) for proteins in unique_prts]
        gene_list = list(set([item for sublist in gene_names for item in sublist]))

        return gene_list

    def get_gene_accessions(self, gene_list: list, species: str = "human"):
        """
        :params gene_list
        """
        mg = mygene.MyGeneInfo()
        gene_accessions = mg.querymany(gene_list, scopes="symbol", species=species, fields="accession")
        gene_accessions_maps = defaultdict(list)
        for obj in gene_accessions:
            if "accession" in obj:
                gene_accessions_maps[obj["query"]].append(obj["accession"])
        return gene_accessions_maps
