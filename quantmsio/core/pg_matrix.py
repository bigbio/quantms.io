from quantmsio.operate.query import Query
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import copy
from quantmsio.core.format import PG_MATRIX

PG_MATRIX_SCHEMA = pa.schema(
    PG_MATRIX,
    metadata={"description": "pg matrix in quantms.io format"},
)

class PgMatrix:
    def __init__(self, feature_path):
        self.db = Query(feature_path)
        
    def get_pep_count(self,df):
        def update_dict(row):
            for p in row['pg_accessions']:
                key = (p, row['reference_file_name'])
                if key in all_dict['peps']:
                    all_dict['peps'][key] += 1
                    all_dict['intensity'][key] += row['intensity']
                else:
                    all_dict['peps'][key] = 1
                    all_dict['intensity'][key] = row['intensity']
        unique_df = df[df['unique']==1][['pg','reference_file_name','intensity']].copy()
        unique_df["peps"] = unique_df.groupby(['pg','reference_file_name']).transform("size")
        unique_df["intensity"] = unique_df.groupby(['pg','reference_file_name'])['intensity'].transform("sum")
        unique_df.set_index(["pg", "reference_file_name"],inplace=True)
        unique_dict = unique_df.to_dict()
        all_dict = copy.deepcopy(unique_dict)
        repeat_df = df[df['unique']==0][['pg_accessions','reference_file_name','intensity']].copy()
        repeat_df.apply(update_dict,axis=1)
        
        return unique_dict,all_dict
    
    def generate_pg_matrix(self):
        def map_count(row):
            res = []
            for name in ['peptides(unique)','peptides(all)']:
                struct = {}
                struct['name'] = name
                if 'unique' in name:
                    struct['value'] = gen_count(row,unique_dict)
                elif 'all' in name:
                    struct['value'] = gen_count(row,all_dict)
                res.append(struct)
            return res
            
        def gen_count(row,map_dict):
            value = []
            ref = row['reference_file_name']
            for p in row['pg_accessions']:
                key = (p,ref)
                if key in map_dict['peps']:
                    value.append(str(map_dict['peps'][key]))
                else:
                    value.append('0')
            return ",".join(value)

        def map_intensity(row):
            res = []
            p = row['pg_accessions']
            ref = row['reference_file_name']
            for name in ['intensity(unique)','intensity(all)']:
                struct = {}
                struct['name'] = name
                key = (p[0],ref)
                if 'unique' in name:
                    if key in unique_dict['intensity']:
                        struct['value'] = unique_dict['intensity'][key]
                    else:
                        struct['value'] = 0
                elif 'all' in name:
                    if key in all_dict['intensity']:
                        struct['value'] = all_dict['intensity'][key]
                    else:
                        struct['value'] = 0
                res.append(struct)
            return res
        
        for refs, df in self.db.iter_file(columns=['pg_accessions','gg_names','unique','quantmsio_version','reference_file_name','intensity']):
            df.loc[:,'pg'] = df['pg_accessions'].str.join(';')
            unique_dict,all_dict = self.get_pep_count(df)
            df.drop_duplicates(subset=['pg'],inplace=True)
            df.loc[:,'peptides'] = df[['pg_accessions','reference_file_name']].apply(lambda row:map_count(row),axis=1)
            df.loc[:,'intensities'] = df[['pg_accessions','reference_file_name']].apply(lambda row: map_intensity(row),axis=1)
            df.loc[:,'first_protein_description'] = None
            yield pa.Table.from_pandas(df, schema=PG_MATRIX_SCHEMA)

    def write_to_file(
        self,
        output_path,
    ):
        pqwriter = None
        for p in self.generate_pg_matrix():
            if not pqwriter:
                pqwriter = pq.ParquetWriter(output_path, p.schema)
            pqwriter.write_table(p)
        if pqwriter:
            pqwriter.close()