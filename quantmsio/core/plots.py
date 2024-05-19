import math
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_distribution_of_ibaq(ibaq_path: str, save_path: str = None, selected_column: str = None) -> None:
    """
    This function plots the distribution of the protein IBAQ values.
    :param ibaq_path: ibaq file path
    :param save_path: save path
    :param selected_column: selected column
    """
    df = pd.read_csv(ibaq_path, sep=None, comment="#", engine="python")
    plt.figure(dpi=500, figsize=(12, 8))
    columns = df.columns
    if selected_column is None:
        if "ribaq" in columns:
            selected_column = "ribaq"
        elif "IbaqLog" in columns:
            selected_column = "IbaqLog"
    if selected_column is not None:
        fig = sns.histplot(data=df[selected_column], stat="frequency", kde=True, color="#209D73")
        fig.set(xlabel=selected_column, ylabel="Frequency")
        fig.set_title("Distribution of IBAQ values using {}".format(selected_column))
        sns.despine(ax=fig, top=True, right=True)
        if save_path:
            fig.figure.savefig(save_path, dpi=500)
        return fig
    else:
        raise ValueError("No IBAQ column found in the ibaq file")


def plot_peptides_of_lfq_condition(psm_parquet_path: str, sdrf_path: str, save_path: str = None) -> None:
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

    df = pd.read_parquet(psm_parquet_path, columns=["reference_file_name"])
    sdrf = pd.read_csv(sdrf_path, sep="\t")
    use_cols = [col for col in sdrf.columns if col.startswith("factor value")]
    use_cols.append("comment[data file]")
    sdrf = sdrf[use_cols]
    sdrf["comment[data file]"] = sdrf["comment[data file]"].apply(lambda x: x.split(".")[0])
    sdrf.rename(columns={"comment[data file]": "reference_file_name"}, inplace=True)
    df = df.merge(sdrf, on="reference_file_name", how="left")
    df.columns = ["reference", "condition"]
    df = df[["condition"]]
    f_count = df["condition"].value_counts()
    f_count.sort_values(ascending=False)
    if len(f_count) < 20:
        i = math.ceil(len(f_count) / 5)
        plt.figure(dpi=500, igsize=(6 * i, 4 * i))
        img = sns.barplot(
            y=f_count.values,
            x=f_count.index,
            hue=f_count.index.astype(str),
            palette="bone_r",
            legend=True,
        )
        img.set(xlabel=None)
        for tick in img.get_xticklabels():
            tick.set_rotation(30)
        sns.despine(ax=img, top=True, right=True)
        if save_path:
            img.figure.savefig(save_path, dpi=500)
        return img
    else:
        df = pd.DataFrame([list(f_count.values)], columns=f_count.index)
        num_subplots = math.ceil(len(f_count) / 20)
        columns_per_subplot = 20
        fig, axes = plt.subplots(nrows=num_subplots, ncols=1, figsize=(12, 4 * num_subplots))
        for i in range(num_subplots):
            start_col = i * columns_per_subplot
            end_col = (i + 1) * columns_per_subplot
            subset_data = df.iloc[:, start_col:end_col]

            sns.barplot(data=subset_data, ax=axes[i])
            axes[i].set_title("Condition vs Number of Peptides {}-{}".format(start_col + 1, end_col))
            axes[i].set(xlabel=None)
            for tick in axes[i].get_xticklabels():
                tick.set_rotation(30)
            sns.despine(ax=axes[i], top=True, right=True)
        plt.tight_layout()
        if save_path:
            fig.figure.savefig(save_path, dpi=500)
        return fig


def plot_intensity_distribution_of_samples(feature_path: str, save_path: str = None, num_samples: int = 10) -> None:
    """
    Plot the distribution of intensities for different samples.
    :param feature_path: path to the feature file
    :param save_path: save path[xxx.svg]
    :param num_samples: number of samples to plot
    """
    df = pd.read_parquet(feature_path, columns=["sample_accession", "intensity"])
    sample_accessions = df["sample_accession"].unique().tolist()
    random.shuffle(sample_accessions)
    if len(sample_accessions) > num_samples:
        sample_accessions = sample_accessions[:num_samples]

    df = df[df["sample_accession"].isin(sample_accessions)]
    df = df[df["intensity"] > 0]
    df["intensity"] = np.log(df["intensity"])

    plt.figure(dpi=500, figsize=(12, 8))
    fig = sns.kdeplot(data=df, x="intensity", hue="sample_accession", palette="Paired", linewidth=2)
    sns.despine(ax=fig, top=True, right=True)
    fig.set(ylabel=None)
    if save_path:
        fig.figure.savefig(save_path, dpi=500)
    return fig


def plot_peptide_distribution_of_protein(feature_path: str, save_path: str = None, num_samples: int = 20) -> None:
    """
    Bar graphs of peptide counts for different samples.
    :param feature_path: path to the feature file
    :param save_path: save path[xxx.svg]
    :param num_samples: number of samples to plot
    """
    df = pd.read_parquet(feature_path, columns=["sample_accession", "sequence"])
    df = df.groupby(["sample_accession"]).agg({"sequence": "count"}).reset_index()
    df.columns = ["sample", "peptides"]
    df = df.sample(frac=1).reset_index(drop=True)
    if len(df) > num_samples:
        df = df.iloc[:20, :]

    max_p = max(df["peptides"]) + 1 / 3 * max(df["peptides"])
    plt.figure(dpi=500, figsize=(12, 8))
    fig = sns.barplot(data=df, x="sample", y="peptides", color="#82C3A3", legend=False, width=0.7)
    fig.set_ylim(0, max_p)
    fig.set(xlabel=None)
    sns.despine(ax=fig, top=True, right=True)
    for tick in fig.get_xticklabels():
        tick.set_rotation(65)
    if save_path:
        fig.figure.savefig(save_path, dpi=500)
    return fig


def plot_intensity_box_of_samples(feature_path: str, save_path: str = None, num_samples: int = 10) -> None:
    """
    Boxplot of peptide intensity distribution for different samples.
    :param feature_path: path to the feature file
    :param save_path: save path[xxx.svg]
    :param num_samples: number of samples to plot
    """
    df = pd.read_parquet(feature_path, columns=["sample_accession", "intensity"])
    sample_accessions = df["sample_accession"].unique().tolist()
    random.shuffle(sample_accessions)

    if len(sample_accessions) > num_samples:
        sample_accessions = sample_accessions[:num_samples]

    df = df[df["sample_accession"].isin(sample_accessions)]
    df = df[df["intensity"] > 0]
    df["intensity"] = np.log(df["intensity"])

    plt.figure(figsize=(12, 8), dpi=500)
    fig = sns.boxplot(
        x="sample_accession",
        y="intensity",
        data=df,
        boxprops=dict(alpha=0.2),
        palette="muted",
        hue="sample_accession",
        fliersize=3,
    )
    fig.set(xlabel=None)
    sns.despine(ax=fig, top=True, right=True)
    for tick in fig.get_xticklabels():
        tick.set_rotation(30)
    if save_path:
        fig.figure.savefig(save_path, dpi=500)
    return fig
