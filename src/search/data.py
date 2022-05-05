"""
Loaders and plotting functions for analysis of sequence data
"""

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import PROTEOMEDATA

# some useful globals
taxcols = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_proteomedata(db):
    """
    Find proteome data and load into a pandas DataFrame.

    :param db: dataset, uniprot or ncbi
    :return: pd.DataFrame
    """
    data = []
    for dataset in PROTEOMEDATA.glob(f"*{db}*.csv"):
        data.append(pd.read_csv(dataset))

    return pd.concat(data)


def isocharge_draw(data, level, name, colorby, n=None):
    """
    Draw fPos-fNeg plots with proteome and protein data. Diagonals represent
    coordinates of identical charge (isocharge lines).

    :param data: a pd.Dataframe containing at least the following columns:
        fPos_x, fNeg_x -> protein features;
        fPos_y,fNeg_y -> proteome features;
        kingdom, phylum, class, order, family, genus -> taxonomic categories
    :param level: one of kingdom, phylum, class, order, family, genus as str
    :param name: a valid taxon name at the level
    :param n: an integer or None limiting the plotting groups
    """
    if n:
        df = data[data[level] == name][:n].copy()
    else:
        df = data[data[level] == name].copy()

    groups = df[colorby].unique()
    palette = sns.color_palette("tab10")
    colors = {g: c for g, c in zip(groups, palette)}
    df["color"] = df[colorby].map(colors)


    fig = plt.figure(figsize=(6,6))
    ax  = fig.subplots(1,1)
    for _, series in df.iterrows():
        ax.plot(
                [series.fPos_x, series.fPos_y],
                [series.fNeg_x, series.fNeg_y],
                linewidth = 0.5,
                marker = "o",
                zorder = 1,
                markersize= 4,
                markerfacecolor = series.color,
                color = series.color,
                markeredgecolor="k"
            )
        ax.scatter(
                    series.fPos_x,
                    series.fNeg_x,
                    marker = "X",
                    zorder = 2,
                    facecolor = series.color,
                    edgecolor = "k",
                    s= 50,
                    label = series.species,

                )
    for i in np.arange(0, 0.25, 0.025):
        ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)

    ax.set_xlabel("fraction positve (KR)")
    ax.set_ylabel("fraction negative (DE)")
    ax.set_title(f"taxonomy({level}) = {name}")
    # ax.legend()
    return fig


def taxonomy_draw(df, x, y, category, bg=True, title=""):
    raise NotImplementedError()
    # fig = plt.figure(figsize=(6,6))
    # ax  = fig.subplots(1,1)
    # # background
    # if bg:
    #     ax.scatter(dfa[x], dfa[y], color = "lightgray", zorder = 0)
    #
    # # foreground
    # for name, subset in df.groupby(category, dropna=False):
    #     name = str(name)
    #     ax.scatter( subset[x], subset[y], label = name,
    #                 edgecolor = "w", marker = "o", linewidth = 0.2, zorder = 2)
    # ax.legend(bbox_to_anchor=(1,1))
    # ax.set_xlabel(x)
    # ax.set_ylabel(y)
    # ax.set_title(title)
    # return fig
