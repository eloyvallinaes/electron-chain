#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from analysis.data import load_carrier_data

# %% codecell
df = load_carrier_data("iter2")


# %% codecell
def draw(data, level, name, n=None):
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

# %% codecell
fig = draw(df, "kingdom", "Archaea")
# fig.savefig("png/relative_Archaea.png", dpi = 300, bbox_inches="tight")

# %% codecell
fig = draw(df, "class", "Halobacteria")
# fig.savefig("png/relative_Halobacteria.png", dpi = 300, bbox_inches="tight")
