import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math


def plot_distribution_of_ibaq(ibaq_path: str, save_path: str, selected_column: str = None) -> None:
    """
    This function plots the distribution of the protein IBAQ values.
    :param ibaq_path: ibaq file path
    :param save_path: save path
    :param selected_column: selected column
    """
    df = pd.read_csv(ibaq_path, sep=None, comment='#', engine='python')
    plt.figure(dpi=500, figsize=(12,8))
    columns = df.columns
    if selected_column is None:
        if 'ribaq' in columns:
            selected_column = 'ribaq'
        elif 'IbaqLog' in columns:
            selected_column = 'IbaqLog'
    if selected_column is not None:
        fig = sns.histplot(data=df[selected_column], stat='frequency', kde=True,color='#209D73')
        fig.set(xlabel=selected_column, ylabel='Frequency')
        fig.set_title('Distribution of IBAQ values using {}'.format(selected_column))
        sns.despine(ax=fig, top=True, right=True)
        fig.figure.savefig(save_path, dpi=500)
    else:
        raise ValueError("No IBAQ column found in the ibaq file")

def plot_peptides_of_lfq_condition(psm_parquet_path: str, sdrf_path: str, save_path:str) -> None:
    """
    this function plots the number of peptides for each condition in a LFQ (Label-Free Quantification) experiment.

    Example Usage
    plot_peptides_of_lfq_condition("psm.parquet", "sdrf.txt", "output.png")
    The function takes three inputs: the path to the PSM (Peptide-Spectrum Match) parquet file, the path to the SDRF
    (Sample and Data Relationship Format) file, and the path to save the output plot.
    It then generates a bar plot showing the number of peptides for each condition in the LFQ
    experiment and saves it as a PNG file.
    :param psm_parquet_path: psm parquet path in lfq
    :param sdrf_path: sdrf path
    :param save_path: save path
    """

    df = pd.read_parquet(psm_parquet_path,columns=["reference_file_name"])
    sdrf = pd.read_csv(sdrf_path,sep='\t')
    use_cols = [col for col in sdrf.columns if col.startswith('factor value')]
    use_cols.append('comment[data file]')
    sdrf = sdrf[use_cols]
    sdrf['comment[data file]'] = sdrf['comment[data file]'].apply(lambda x: x.split('.')[0])
    sdrf.rename(columns={'comment[data file]': "reference_file_name"}, inplace=True)
    df = df.merge(sdrf,on="reference_file_name",how='left')
    df.columns = ['reference','condition']
    df = df[['condition']]
    f_count = df['condition'].value_counts()
    f_count.sort_values(ascending=False)
    if len(f_count) < 20:
        i = math.ceil(len(f_count)/5)
        plt.figure(dpi=500, igsize=(6*i,4*i))
        img = sns.barplot(y=f_count.values, x=f_count.index, hue=f_count.index.astype(str),
                          palette="bone_r", legend=True)
        img.set(xlabel=None)
        for tick in img.get_xticklabels():
            tick.set_rotation(30)
        sns.despine(ax=img, top=True, right=True)
        img.figure.savefig(save_path,dpi=500)
    else:
        df = pd.DataFrame([list(f_count.values)], columns=f_count.index)
        num_subplots = math.ceil(len(f_count)/20)
        columns_per_subplot = 20
        fig, axes = plt.subplots(nrows=num_subplots, ncols=1, figsize=(12,4*num_subplots))
        for i in range(num_subplots):
            start_col = i * columns_per_subplot
            end_col = (i + 1) * columns_per_subplot
            subset_data = df.iloc[:, start_col:end_col]

            sns.barplot(data=subset_data, ax=axes[i])
            axes[i].set_title("Condition vs Number of Peptides {}-{}".format(start_col+1, end_col))
            axes[i].set(xlabel=None)
            for tick in axes[i].get_xticklabels():
                tick.set_rotation(30)
            sns.despine(ax=axes[i], top=True, right=True)
        plt.tight_layout()
        fig.figure.savefig(save_path, dpi=500)
