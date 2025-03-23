from pathlib import Path
from typing import Union, Optional

import pandas as pd
import anndata as ad


def pivot_wider(
    df: pd.DataFrame,
    row_name: str,
    col_name: str,
    values: str,
    fillna=False,
) -> pd.DataFrame:
    """
    Create a matrix from a DataFrame given the row, column, and value columns.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame in long format.
    row_name : str
        The column name to use as row labels (e.g., sample_ids).
    col_name : str
        The column name to use as column labels (e.g., protein_names).
    values : str
        The column name to use as cell values (e.g., expression_values).
    fillna : Optional[Union[bool, int, float]]
        Value to fill NaN. If True, fill NaN with 0. If False or None, leave NaN as is.
        If a number is provided, use that value.

    Returns
    -------
    pd.DataFrame
        A pivot table (matrix) with specified rows, columns, and values.

    Examples
    --------
    >>> df_matrix =  pivot_wider(combined_df,
                    row_name='SampleID',
                    col_name='ProteinName',
                    values='Ibaq',
                    fillna=False)
    """
    # Check if the provided columns exist in the DataFrame
    missing_columns = {row_name, col_name, values} - set(df.columns)
    if missing_columns:
        raise ValueError(f"Columns {missing_columns} not found in the DataFrame.")

    # Check for duplicate combinations
    duplicates = df.groupby([row_name, col_name]).size()
    if (duplicates > 1).any():
        raise ValueError(
            f"Found duplicate combinations of {row_name} and {col_name}. "
            "Use an aggregation function to handle duplicates."
        )

    # Use pivot_table to create the matrix
    matrix = df.pivot_table(
        index=row_name, columns=col_name, values=values, aggfunc="first"
    )

    # Simplified NaN handling
    if fillna is True:  # Fill with 0 if True
        matrix = matrix.fillna(0)
    elif fillna not in [None, False]:  # Fill if a specific value is provided
        matrix = matrix.fillna(fillna)

    return matrix


def get_condition_map(df: pd.DataFrame) -> dict:
    m2 = df.drop_duplicates(subset=["sample_accession"])
    m2.set_index("sample_accession", inplace=True)
    s2 = m2.to_dict()["condition"]
    return s2


class Combiner:
    def __init__(self):
        self.combined_adata: Optional[ad.AnnData] = None

    @staticmethod
    def transform_to_adata(df: pd.DataFrame) -> ad.AnnData:
        express_matrix = pivot_wider(
            df, row_name="sample_accession", col_name="protein", values="ibaq"
        )
        express2_matrix = pivot_wider(
            df,
            row_name="sample_accession",
            col_name="protein",
            values="ibaq_normalized",
        )
        adata = ad.AnnData(
            X=express_matrix.to_numpy(),
            obs=express_matrix.index.to_frame(),
            var=express_matrix.columns.to_frame(),
        )
        condition_map = get_condition_map(df)
        adata.obs["condition"] = adata.obs["sample_accession"].map(condition_map)
        adata.layers["ibaq_normalized"] = express2_matrix.to_numpy()
        return adata

    def combine_adata(self, adata: ad.AnnData, axis=0, join="outer") -> None:
        if self.combined_adata is None:
            self.combined_adata = adata
        else:
            self.combined_adata = ad.concat(
                [self.combined_adata, adata], axis=axis, join=join
            )
        return self.combined_adata

    def save_adata(self, output_path: Union[Path, str]) -> None:
        if self.combined_adata is not None:
            self.combined_adata.write(output_path)

    def filter_contion(self, condition: str) -> Optional[ad.AnnData]:
        if self.combined_adata is None:
            return None
        else:
            return self.combined_adata[self.combined_adata.obs.condition == condition]

    def to_df(self, layer=None) -> Optional[pd.DataFrame]:
        if self.combined_adata is None:
            return None
        else:
            return self.combined_adata.to_df(layer=layer)
