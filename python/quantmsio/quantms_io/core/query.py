import os
import re
import duckdb


def check_string(re_exp, strings):
    res = re.search(re_exp, strings)
    if res:
        return True
    else:
        return False


class FeatureQuery:

    def __init__(self, parquet_path: str):
        if os.path.exists(parquet_path):
            self.parquet_db = duckdb.connect()
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path))
        else:
            raise FileNotFoundError(f'the file {parquet_path} does not exist.')

    def get_unique_peptides(self):
        """
        return: A list of deduplicated peptides.
        """
        feature_db = self.feature_db
        unique_peps = self.parquet_db.sql(f"SELECT DISTINCT sequence FROM parquet_db").df()

        return unique_peps['sequence'].tolist()

    def get_unique_proteins(self):
        """
        return: A list of deduplicated proteins.
        """

        unique_prts = self.parquet_db.sql(f"SELECT DISTINCT protein_accessions FROM parquet_db").df()

        return unique_prts['protein_accessions'].tolist()

    def query_peptide(self, peptide: str):
        """
        peptide: Peptide that need to be queried.
        return: A DataFrame of all information about query peptide.
        """

        if check_string('^[A-Z]+$', peptide):
            return self.parquet_db.sql(f"SELECT * FROM parquet_db WHERE sequence ='{peptide}'").df()
        else:
            return KeyError('Illegal peptide!')


    def query_protein(self, protein: str):
        """
        protein: Protein that need to be queried.
        return: A DataFrame of all information about query protein.
        """
        if check_string('^[A-Z]+', protein):
            return self.parquet_db.sql(f"SELECT * FROM parquet_db WHERE protein_accessions ='{protein}'").df()
        else:
            return KeyError('Illegal protein!')