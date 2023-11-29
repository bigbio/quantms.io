import click
import pandas as pd
import pyarrow.parquet as pq


def read_large_parquet(parquet_path: str, batch_size: int = 100000):
    """_summary_

    :param parquet_path: _description_
    :param batch_size: _description_, defaults to 100000
    :yield: _description_
    """
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batch_df = batch.to_pandas()
        yield batch_df


@click.command()
@click.option("--feature_parquet", help="Parquet file with features", required=True)
@click.option("--output_file", help="Output file", required=True)
@click.option("--accession", help="Accession", required=True)
@click.option("--batch_size", help="Batch size", default=100000)
def extract_protein_qvalues(feature_parquet: str, output_file: str, accession: str, batch_size: int = 100000):

    chunks = read_large_parquet(feature_parquet, batch_size=batch_size)
    q_df = pd.DataFrame()
    for feature_df in chunks:
        feature_df = feature_df[['protein_accessions', 'protein_global_qvalue', 'is_decoy', 'condition']]
        feature_df = feature_df.explode('protein_accessions')
        q_df = pd.concat([q_df, feature_df])
        q_df.drop_duplicates(inplace=True)
        q_df['dataset_accession'] = accession
    print(q_df)
    q_df.to_csv(output_file, index=False, mode="w")

if __name__ == '__main__':
    try:
        extract_protein_qvalues()
    except Exception as e:
        print(e)
