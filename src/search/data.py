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
    Find proteome data and load into a pandas DataFrame. For ncbi data, create
    Name column for consistent merging with targets.

    :param db: dataset, <uniprot|ncbi>
    :return: pd.DataFrame
    """
    data = []
    for dataset in PROTEOMEDATA.glob(f"*{db}*.csv"):
        data.append(pd.read_csv(dataset))

    df = pd.concat(data)
    if db == "ncbi":
        df["Name"] = df.AssemblyID

    return df


def isocharge_artist(ax, df, label = None, color = "tab:blue", diagonals=True):
    for _, series in df.iterrows():
        ax.plot(
                [series.fPos_x, series.fPos_y],
                [series.fNeg_x, series.fNeg_y],
                linewidth = 0.5,
                marker = "",
                zorder = 1,
                color = color,
            )
    # circles
    ax.scatter(
                df.fPos_y, # proteome coordinates
                df.fNeg_y,
                marker = "o",
                zorder = 1,
                facecolor = color,
                edgecolor = "k",
                linewidth = 0.2,
                s= 50,

            )
    # crosses
    ax.scatter(
                df.fPos_x, # protein coordinates
                df.fNeg_x,
                marker = "X",
                zorder = 2,
                facecolor = color,
                edgecolor = "k",
                linewidth = 0.2,
                s= 50,
                label = label,
            )

    if diagonals:
        # iso-charge diagonals
        for i in np.arange(0, 0.25, 0.025):
            ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)

    return ax


def target_proteome_merge(targets, proteomes):
    """
    Combine protein-proteome dataframes.
    """
    return targets.merge(
                    proteomes.drop(taxcols, axis = 1),
                    left_on="proteomecode",
                    right_on="Name",
                    how = "left"
                )


def isocharge_draw(ax, df, rank = "", name = "",
                   lower = "", groups = 10, colors = None):
    """
    Draw fPos-fNeg plots with proteome and protein data. Diagonals represent
    coordinates of identical charge (isocharge lines).

    :param data: a pd.Dataframe containing at least the following columns:
        fPos_x, fNeg_x -> protein features;
        fPos_y,fNeg_y -> proteome features;
        kingdom, phylum, class, order, family, genus -> taxonomic categories
    :rank param: str, one of <kingdom|phylum|class|order|family|genus>
    :name param: str, the group at rank to split
    :lower param: str, lower rank or any categorical to group and color data
    :groups param: int or list, if int, the number of groups to consider at lower rank,
                   if list, the names of the groups
    :background param: bool, if True, include all observation in gray
    :colors param: a dictionary mapping groups to colors, if None autogenerate
    """

    # groups
    if type(groups) == int:
        if rank == "" or name == "":
            data = df
        else:
            data = df[df[rank] == name]

        if lower == "":
            groups = []
        else:
            groups = data[lower].value_counts().index[:groups]

    # color for each group
    if not colors:
        colors = make_cmap(groups)

    if len(groups) == 0:
        isocharge_artist(ax, data)
    else:
        for g in groups:
            d = data[data[lower] == g]
            isocharge_artist(ax, d, label = g, color = colors[g])

    ax.set_xlabel("fraction positve (KR)")
    ax.set_ylabel("fraction negative (DE)")
    ax.set_title(f"taxonomy({rank}) = {name}")
    return ax


def taxonomy_draw( ax, df, x, y, rank = "kingdom", name = "Bacteria",
                   lower = "phylum", groups = 10, background = True,
                   colors = None, **kwargs
            ):
    """
    Draw 2D scatter and color by taxonomic group.

    args
    :ax param: a plt.Axes instance where to plot
    :df param: a pd.DataFrame instance with columns x,y and taxonomical
    :rank param: str, one of <kingdom|phylum|class|order|family|genus>
    :name param: str, the group at rank to split
    :lower param: str, lower rank or any categorical to group and color data
    :groups param: int or list, if int, the number of groups to consider at lower rank,
                   if list, the names of the groups
    :background param: bool, if True, include all observation in gray
    :colors param: a dictionary mapping groups to colors, if None autogenerate

    kwargs
    pass to ax.scatter to modify marker properties eg. s (size)
    """
    # background
    if background:
        ax.scatter(df[x], df[y], c="lightgray", s=20, zorder=1, label="")

    # groups
    if type(groups) == int:
        subset = df[df[rank] == name]
        groups = subset[lower].value_counts().index[:groups]

    # color for each group
    if not colors:
        colors = make_cmap(groups)

    for group in groups:
        data = df[df[lower] == group]
        ax.scatter(
                    x=data[x], y=data[y], c = colors[group], label = group,
                    zorder = 10, edgecolor = "w", linewidth = 0.3, **kwargs
                )
    ax.set_title(f"{rank} = {name}; into {lower}")
    return ax


def make_cmap(groups):
    """
    Associate tab10 colors to names in a list of groups.
    """
    palette = sns.color_palette('tab10', len(groups))
    colordict = {x: palette[i] for i, x in enumerate(groups)}
    return colordict
