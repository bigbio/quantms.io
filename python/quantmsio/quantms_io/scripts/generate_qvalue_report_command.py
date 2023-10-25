import click
import pandas as pd


def read_large_parquet(parquet_path: str, batch_size: int = 100000):
    """_summary_

    :param parquet_path: _description_
    :param batch_size: _description_, defaults to 100000
    :yield: _description_
    """
    reg_df = pd.read_parquet(parquet_path, engine='pyarrow', columns=['protein_accessions', 'protein_global_qvalue', 'is_decoy', 'condition'])
    return reg_df


@click.command()
@click.option("--feature_parquet", help="Parquet file with features", required=True)
@click.option("--output_file", help="Output file", required=True)
def extract_protein_qvalues(feature_parquet: str, output_file: str):

    df_qvalues = read_large_parquet(feature_parquet)
    df_qvalues.drop_duplicates(inplace=True)
    df_qvalues.to_csv(output_file, index=False, mode='w')

if __name__ == '__main__':
    try:
        extract_protein_qvalues()
    except Exception as e:
        print(e)
