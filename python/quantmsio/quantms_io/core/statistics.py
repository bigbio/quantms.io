import os
import re
from abc import ABC

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
import random


class Statistics(ABC):

    def get_number_of_proteins(self):
        pass

    def get_number_of_peptides(self):
        pass

    def get_number_of_samples(self):
        pass


class IbaqStatistics(Statistics):

    def __init__(self, ibaq_path: str) -> None:
        self.ibaq_path = ibaq_path
        self.ibaq_db = pd.read_csv(ibaq_path, sep=None, comment='#', engine='python')

    def get_number_of_proteins(self) -> int:
        if 'ProteinName' in self.ibaq_db.columns:
            return len(self.ibaq_db['ProteinName'].unique())
        elif 'protein' in self.ibaq_db.columns:
            return len(self.ibaq_db['protein'].unique())
        else:
            raise ValueError("No protein column found in the ibaq file")

    def get_number_of_samples(self) -> int:
        if 'SampleID' in self.ibaq_db.columns:
            return len(self.ibaq_db['SampleID'].unique())
        else:
            raise ValueError("No SampleID column found in the ibaq file")


class ParquetStatistics(Statistics):

    def __init__(self, parquet_path: str) -> None:
        if os.path.exists(parquet_path):
            self.parquet_db = duckdb.connect()
            self.parquet_db = self.parquet_db.execute(
                "CREATE VIEW parquet_db AS SELECT * FROM parquet_scan('{}')".format(parquet_path))
        else:
            raise FileNotFoundError(f'the file {parquet_path} does not exist.')

    def get_number_of_peptides(self) -> int:
        count = self.parquet_db.sql(f"SELECT COUNT(DISTINCT sequence) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_peptidoforms(self) -> int:
        count = self.parquet_db.sql(f"SELECT COUNT(DISTINCT peptidoform) FROM parquet_db").fetchone()[0]
        return count

    def get_number_of_samples(self) -> int:
        count = self.parquet_db.sql(f"SELECT COUNT(DISTINCT sample_accession) FROM parquet_db").fetchone()[0]
        return count

    def get_number_proteins(self) -> int:
        """
        This method reads the protein accessions from a parquet file and return the number of unique accessions.
        This method is not accurate. It needs to be refined.
        :return: number of unique proteins
        """
        protein_ids = self.parquet_db.sql(f"SELECT DISTINCT protein_accessions FROM parquet_db").fetchall()
        # This probalby needs to be refined.
        protein_ids = [item for sublist in protein_ids for item in sublist]
        protein_ids = [item for sublist in protein_ids for item in sublist]
        protein_ids = set(protein_ids)
        return len(protein_ids)



def check_string(re_exp, strings):
    res = re.search(re_exp, strings)
    if res:
        return True
    else:
        return False


class Statistic:
    __slot__ = ['feature_db', 'ibaq_db']

    def __init__(self, feature_path: str = None, ibaq_path: str = None):

        self.feature_db = self.load_db(feature_path)
        self.ibaq_db = self.load_db(ibaq_path)

    def load_db(self, path: str):
        """
        Load the database.
        """
        if path is not None:
            if os.path.exists(path):
                if path.endswith(".parquet"):
                    return duckdb.read_parquet(path)
                else:
                    skip_rows = get_skip_rows(path)
                    return duckdb.read_csv(path, sep='\t', skiprows=skip_rows)
            else:
                return FileNotFoundError(f'the file {path} does not exist.')
        else:
            return None

    def get_unique_peptides(self):
        """
        return: A list of deduplicated peptides.
        """
        feature_db = self.feature_db
        unique_peps = duckdb.sql(f"SELECT DISTINCT sequence FROM feature_db").df()

        return unique_peps['sequence'].tolist()

    def get_unique_proteins(self):
        """
        return: A list of deduplicated proteins.
        """
        ibaq_db = self.ibaq_db
        unique_prts = duckdb.sql(f"SELECT DISTINCT protein FROM ibaq_db").df()

        return unique_prts['protein'].tolist()

    def query_peptide(self, peptide: str):
        """
        peptide: Peptide that need to be queried.
        return: A DataFrame of all information about query peptide.
        """

        feature_db = self.feature_db
        if check_string('^[A-Z]+$', peptide):
            return duckdb.sql(f"SELECT * FROM feature_db WHERE sequence ='{peptide}'").df()
        else:
            return KeyError('Illegal peptide!')

    def query_protein(self, protein: str):
        """
        protein: Protein that need to be queried.
        return: A DataFrame of all information about query protein.
        """

        ibaq_db = self.ibaq_db
        if check_string('^[A-Z]+', protein):
            return duckdb.sql(f"SELECT * FROM ibaq_db WHERE protein ='{protein}'").df()
        else:
            return KeyError('Illegal protein!')

    def plot_peptide_distribution_of_protein(self):
        """
        Bar graphs of peptide counts for different samples.
        """
        feature_db = self.feature_db
        df = duckdb.sql(f"SELECT sample_accession, COUNT(sequence) FROM feature_db GROUP BY sample_accession").df()
        df.columns = ['sample', 'peptides']
        df = df.sample(frac=1).reset_index(drop=True)
        if len(df) > 20:
            return df.iloc[:20, :].plot.bar(figsize=(12, 8), x='sample', y='peptides',
                                            title='number of peptides for different samples', width=0.7,
                                            color='#82C3A3', rot=65, ylim=(
                0, max(df.loc[:20, 'peptides']) + 1 / 2 * max(df.loc[:20, 'peptides'])))
        else:
            return df.plot.bar(figsize=(12, 8), x='sample', y='peptides',
                               title='number of peptides for different samples', color='#82C3A3', rot=65)

    def plot_intensity_distribution_of_samples(self):
        """
        Kde of peptide intensity distribution for different samples.
        """
        feature_db = self.feature_db
        sample_accessions = duckdb.sql(f"SELECT DISTINCT sample_accession FROM feature_db").df()[
            'sample_accession'].tolist()
        random.shuffle(sample_accessions)
        if len(sample_accessions) > 10:
            sample_accessions = sample_accessions[:10]
        df = pd.DataFrame()
        for sample in sample_accessions:
            df_sample = duckdb.sql(f"SELECT intensity FROM feature_db WHERE sample_accession='{sample}'").df()
            df_sample = df_sample[df_sample['intensity'] > 0]
            df_sample['intensity'] = np.log(df_sample['intensity'])
            df_sample.columns = [sample]
            df = pd.concat([df, df_sample], axis=1)

        return df.plot.kde(figsize=(12, 8), linewidth=2, legend=False)


    def plot_intensity_box_of_samples(self):
        """
        Boxplot of peptide intensity distribution for different samples.
        """
        feature_db = self.feature_db
        sample_accessions = duckdb.sql(f"SELECT DISTINCT sample_accession FROM feature_db").df()[
            'sample_accession'].tolist()
        random.shuffle(sample_accessions)
        if len(sample_accessions) > 10:
            sample_accessions = sample_accessions[:10]
        df = pd.DataFrame()
        for sample in sample_accessions:
            df_sample = duckdb.sql(f"SELECT sample_accession,intensity FROM db WHERE sample_accession='{sample}'").df()
            df_sample = df_sample[df_sample['intensity'] > 0]
            df_sample['intensity'] = np.log(df_sample['intensity'])
            df = pd.concat([df, df_sample], axis=0)
        plt.figure(figsize=(12, 8), dpi=500)
        chart = sns.boxplot(
            x='sample_accession',
            y='intensity',
            data=df,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )
        chart.set_xticklabels(labels=sample_accessions, rotation=65)

        return chart
