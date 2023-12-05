import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
import random

def get_skip_rows(path):
    skip_rows = 0
    with open(path) as f:
        line = f.readline()
        while line.startswith("#"):
            skip_rows += 1
            line = f.readline()
        return skip_rows

def check_string(re_exp, strings):
    res = re.search(re_exp, strings)
    if res:
        return True
    else:
        return False
    
def check_exist(label):
    def check_db_exist(func):
        def wrapper(self, *args, **kwargs):
            if label=='feature_db':
                if self.feature_db is None:
                    return LookupError('The database is not loaded.')
            elif label=='ibaq_db':
                if self.load_db is None:
                    return LookupError('The database is not loaded.')
            
            res = func(self, *args, **kwargs)
            return res
        return wrapper
    return check_db_exist
    
class Statistic:
    
    __slot__ = ['feature_db','ibaq_db']
    
    def __init__(self,feature_path:str=None,ibaq_path:str=None):
        
        self.feature_db = self.load_db(feature_path)
        self.ibaq_db = self.load_db(ibaq_path)
        
    def load_db(self,path:str):
        """
        Load the database.
        """
        if path is not None:
            if os.path.exists(path):
                if path.endswith(".parquet"):
                    return duckdb.read_parquet(path)
                else:
                    skip_rows = get_skip_rows(path)
                    return duckdb.read_csv(path,sep='\t',skiprows=skip_rows)
            else:
                return FileNotFoundError(f'the file {path} does not exist.')
        else:
            return None
                
    @check_exist('feature_db')
    def get_unique_peptides(self):
        """
        return: A list of deduplicated peptides.
        """
        feature_db = self.feature_db
        unique_peps = duckdb.sql(f"SELECT DISTINCT sequence FROM feature_db").df()
        
        return unique_peps['sequence'].tolist()

    @check_exist('ibaq_db')
    def get_unique_proteins(self):
        """
        return: A list of deduplicated proteins.
        """
        ibaq_db = self.ibaq_db
        unique_prts = duckdb.sql(f"SELECT DISTINCT protein FROM ibaq_db").df()
        
        return unique_prts['protein'].tolist()
 
    @check_exist('feature_db')
    def query_peptide(self,peptide:str):
        """
        peptide: Peptide that need to be queried.
        return: A DataFrame of all information about query peptide.
        """

        feature_db = self.feature_db
        if check_string('^[A-Z]+$',peptide):
            return duckdb.sql(f"SELECT * FROM feature_db WHERE sequence ='{peptide}'").df()
        else:
            return KeyError('Illegal peptide!')

    @check_exist('ibaq_db')
    def query_protein(self,protein:str):
        """
        protein: Protein that need to be queried.
        return: A DataFrame of all information about query protein.
        """

        ibaq_db = self.ibaq_db
        if check_string('^[A-Z]+',protein):
            return duckdb.sql(f"SELECT * FROM ibaq_db WHERE protein ='{protein}'").df()
        else:
            return KeyError('Illegal protein!')
    
    @check_exist('feature_db')
    def plot_peptide_distribution_of_protein(self):
        """
        Bar graphs of peptide counts for different samples.
        """
        feature_db = self.feature_db
        df = duckdb.sql(f"SELECT sample_accession, COUNT(sequence) FROM feature_db GROUP BY sample_accession").df()
        df.columns = ['sample','peptides']
        df = df.sample(frac=1).reset_index(drop=True)
        if len(df) > 20:
            return df.iloc[:20,:].plot.bar(figsize=(12,8),x='sample',y='peptides',title='number of peptides for different samples',width=0.7,color='#82C3A3',rot=65,ylim=(0,max(df.loc[:20,'peptides'])+1/2*max(df.loc[:20,'peptides'])))
        else:
            return df.plot.bar(figsize=(12,8),x='sample',y='peptides',title='number of peptides for different samples',color='#82C3A3',rot=65)
    
    @check_exist('feature_db')
    def plot_intensty_distribution_of_samples(self):
        """
        Kde of peptide intensity distribution for different samples.
        """
        feature_db = self.feature_db
        sample_accessions = duckdb.sql(f"SELECT DISTINCT sample_accession FROM feature_db").df()['sample_accession'].tolist()
        random.shuffle(sample_accessions)
        if len(sample_accessions) > 10:
            sample_accessions = sample_accessions[:10]
        df = pd.DataFrame()
        for sample in sample_accessions:
            df_sample = duckdb.sql(f"SELECT intensity FROM feature_db WHERE sample_accession='{sample}'").df()
            df_sample = df_sample[df_sample['intensity']>0]
            df_sample['intensity'] = np.log(df_sample['intensity'])
            df_sample.columns = [sample]
            df = pd.concat([df,df_sample],axis=1)
        
        return df.plot.kde(figsize=(12,8),linewidth=2, legend=False)
    
    @check_exist('feature_db')
    def plot_intensty_box_of_samples(self):
        """
        Boxplot of peptide intensity distribution for different samples.
        """
        feature_db = self.feature_db
        sample_accessions = duckdb.sql(f"SELECT DISTINCT sample_accession FROM feature_db").df()['sample_accession'].tolist()
        random.shuffle(sample_accessions)
        if len(sample_accessions) > 10:
            sample_accessions = sample_accessions[:10]
        df = pd.DataFrame()
        for sample in sample_accessions:
            df_sample = duckdb.sql(f"SELECT sample_accession,intensity FROM db WHERE sample_accession='{sample}'").df()
            df_sample = df_sample[df_sample['intensity']>0]
            df_sample['intensity'] = np.log(df_sample['intensity'])
            df = pd.concat([df,df_sample],axis=0)
        plt.figure(figsize=(12,8),dpi=500)
        chart = sns.boxplot(
            x='sample_accession',
            y='intensity',
            data=df,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )
        chart.set_xticklabels(labels =sample_accessions, rotation = 65)
        
        return chart
