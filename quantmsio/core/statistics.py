import os
from abc import ABC

import duckdb
import pandas as pd


class Statistics(ABC):

    def get_number_of_proteins(self):
        pass

    def get_number_of_peptides(self):
        pass

    def get_number_of_samples(self):
        pass

    def get_number_of_peptidoforms(self):
        pass

    def get_number_msruns(self):
        pass


class IbaqStatistics(Statistics):

    def __init__(self, ibaq_path: str) -> None:
        self.ibaq_path = ibaq_path
        self.ibaq_db = pd.read_csv(ibaq_path, sep=None, comment="#", engine="python")

    def get_number_of_proteins(self) -> int:
        if "ProteinName" in self.ibaq_db.columns:
            return len(self.ibaq_db["ProteinName"].unique())
        elif "protein" in self.ibaq_db.columns:
            return len(self.ibaq_db["protein"].unique())
        else:
            raise ValueError("No protein column found in the ibaq file")

    def get_number_of_samples(self) -> int:
        if "SampleID" in self.ibaq_db.columns:
            return len(self.ibaq_db["SampleID"].unique())
        elif "sample_accession" in self.ibaq_db.columns:
            return len(self.ibaq_db["sample_accession"].unique())
        else:
            raise ValueError("No SampleID column found in the ibaq file")


class ParquetStatistics(Statistics):

    def __init__(self, parquet_path: str) -> None:
        if os.path.exists(parquet_path):
            self.parquet_db = duckdb.connect()
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path)
            )
        else:
            raise FileNotFoundError(f"the file {parquet_path} does not exist.")

    def get_number_of_peptides(self) -> int:
        count = self.parquet_db.sql("SELECT COUNT(DISTINCT sequence) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_peptidoforms(self) -> int:
        count = self.parquet_db.sql("SELECT COUNT(DISTINCT peptidoform) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_samples(self) -> int:
        count = self.parquet_db.sql("SELECT COUNT(DISTINCT sample_accession) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_proteins(self) -> int:
        """
        This method reads the protein accessions from a parquet file and return the number of unique accessions.
        This method is not accurate. It needs to be refined.
        :return: number of unique proteins
        """
        protein_ids = self.parquet_db.sql("SELECT DISTINCT protein_accessions FROM parquet_db").fetchall()
        # This probalby needs to be refined.
        protein_ids = [item for sublist in protein_ids for item in sublist[0]]
        protein_ids = set(protein_ids)
        return len(protein_ids)

    def get_number_msruns(self) -> int:
        count = self.parquet_db.sql("SELECT COUNT(DISTINCT reference_file_name) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_psms(self) -> int:
        """
        If the file is a psm file, it will return the number of psms
        :return: numbers of psms
        """
        count = self.parquet_db.sql("SELECT COUNT(*) FROM parquet_db").fetchone()[0]
        return count
