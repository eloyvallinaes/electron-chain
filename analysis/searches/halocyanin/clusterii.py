#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from analysis.data import load_carrier_data, load_ncbi_data

# %% codecell
df = load_carrier_data("halocyanin", "ClusterII", 0, "ncbi")
ncbi = load_ncbi_data()
dfa = ncbi[ncbi.kingdom == "Archaea"]


# %% codecell
def isocharge_draw(data, level, name, n=None):
    if n:
        df = data[data[level] == name][:n]
    else:
        df = data[data[level] == name]

    fig = plt.figure(figsize=(6,6))
    ax  = fig.subplots(1,1)
    for series in df.itertuples():
        ax.plot(
                [series.fPos_x, series.fPos_y],
                [series.fNeg_x, series.fNeg_y],
                label = series.species,
                linewidth = 0.5,
                marker = "o",
                zorder = 1,
                markersize= 4,
                markeredgecolor="k"
            )
        ax.scatter(
                    series.fPos_x,
                    series.fNeg_x,
                    marker = "X",
                    zorder = 2,
                    edgecolor = "k",
                    s= 50
                )
    for i in np.arange(0, 0.25, 0.025):
        ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)
        # ax.plot([0,.2], [.1,.3], linestyle = "--", color= "k")

    ax.set_xlabel("fraction positve (KR)")
    ax.set_ylabel("fraction negative (DE)")
    ax.set_title(f"taxonomy({level}) = {name}")
    ax.set_xlim([0,.15])
    # ax.set_ylim([0,20])


    return fig

def taxonomy_draw(df, x, y, category, bg=True, title=""):
    fig = plt.figure(figsize=(6,6))
    ax  = fig.subplots(1,1)
    # background
    if bg:
        ax.scatter(dfa[x], dfa[y], color = "lightgray", zorder = 0)

    # foreground
    for name, subset in df.groupby(category, dropna=False):
        name = str(name)
        ax.scatter( subset[x], subset[y], label = name,
                    edgecolor = "w", marker = "o", linewidth = 0.2, zorder = 2)
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    return fig


# %% codecell
fig = isocharge_draw(df, "kingdom", "Archaea")
# fig.savefig("png/isocharge_all.png", dpi = 300, bbox_inches="tight")

# %% codecell
fig = isocharge_draw(df, "class", "Candidatus Poseidoniia")
# fig.savefig("png/isocharge_poseidoniia.png", dpi = 300, bbox_inches="tight")

# %% codecell
fig = taxonomy_draw(df, "NCD1000", "MWkDa", "class", bg = True,
                    title="PROTEOME properties")
# fig.savefig("png/mwkda_ncd1000_proteomes.png", dpi = 300, bbox_inches="tight")

# %% codecell
fig = taxonomy_draw(df, "ncd1000", "mwkda", "class", bg = False,
                    title="protein properties")
# fig.savefig("png/mwkda_ncd1000_proteins.png", dpi = 300, bbox_inches="tight")


# %% codecell
# List for Nannan
selection =( df["class"] == "Candidatus Poseidoniia" ) & ( df.mwkda < 40 )
df[selection][['ncd1000', 'mwkda', "MWkDa", "NCD1000", "class", "upacc", "name"]].to_csv("cii_halocyanin_targets.csv", index=False)
