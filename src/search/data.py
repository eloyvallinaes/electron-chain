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
    colors['nan'] = "lightgray"


    fig = plt.figure(figsize=(6,6))
    ax  = fig.subplots(1,1)
    for _, series in df.iterrows():
        ax.plot(
                [series.fPos_x, series.fPos_y],
                [series.fNeg_x, series.fNeg_y],
                linewidth = 0.5,
                marker = "",
                zorder = 1,
                color = "k",
            )
    for g, subset in df.groupby(colorby):
        ax.scatter(
                    subset.fPos_y, # proteome coordinates
                    subset.fNeg_y,
                    marker = "o",
                    zorder = 1,
                    facecolor = colors[g],
                    edgecolor = "k",
                    s= 50,
                    label = g,

                )
        ax.scatter(
                    subset.fPos_x, # protein coordinates
                    subset.fNeg_x,
                    marker = "X",
                    zorder = 2,
                    facecolor = colors[g],
                    edgecolor = "k",
                    s= 50,
                    label = g,

                )
    for i in np.arange(0, 0.25, 0.025):
        ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)

    ax.set_xlabel("fraction positve (KR)")
    ax.set_ylabel("fraction negative (DE)")
    ax.set_title(f"taxonomy({level}) = {name}")
    # ax.legend()
    return fig


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
    :lower param: str, the lower rank into which the data is split
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
    # lightgray RGB code, supress warning about ragged ndarrays
    palette.append((211/254, 211/254, 211/254))
    colordict = {x: [palette[i]] for i, x in enumerate(groups)}
    return colordict
