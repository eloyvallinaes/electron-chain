#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import HALO_DIR, CYTO_DIR
from src.search import filter
from src.search.data import load_proteomedata, isocharge_draw, taxcols, taxonomy_draw, make_cmap

# %% codecell
with open(CYTO_DIR / "bacteria/uniprot/iter1.csv", "r") as cyto:
    reader = csv.DictReader(cyto)
    hits = [h for h in reader]
    hits = filter.top_per_proteome(hits)
    hits["ncd1000"] = hits.ncd * 1000
    hits["mwkda"] = hits.mass / 1000

up = load_proteomedata("uniprot").drop(taxcols, axis = "columns")
df = hits.merge(up, left_on="proteomecode", right_on="Name", how="inner")


# %% codecell
params = dict(rank = "phylum", name = "Proteobacteria", lower="class", groups = 5)
fig = plt.figure(figsize=(6,6))
ax = fig.subplots(1,1)
taxonomy_draw(ax, hits, "ncd1000", "mwkda", **params)
ax.legend()

fig = plt.figure(figsize=(6,6))
ax = fig.subplots(1,1)
taxonomy_draw(ax, hits, "fPos", "fNeg", **params)
ax.legend()


# %% codecell
fig = plt.figure(figsize=(6,12))
axs  = fig.subplots(2,1)
#ax0
taxonomy_draw(axs[0], df, "fPos_y", "fNeg_y", rank = "class", name="Gammaproteobacteria", lower="order", groups=4)
for i in np.arange(0, 0.25, 0.025):
    axs[0].plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)
axs[0].legend()
# ax1
taxonomy_draw(axs[1], df, "fPos_x", "fNeg_x", rank = "class", name="Gammaproteobacteria", lower="order", groups=4)
for i in np.arange(0, 0.25, 0.025):
    axs[1].plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)


# %% codecell
refined = df[df.name.str.fullmatch("[Cc]ytochrome [cC]")]


# %% codecell
fig = plt.figure(figsize=(6,12))
axs  = fig.subplots(2,1)
params = dict(rank = "kingdom", name="Bacteria", lower="phylum", groups=4)
#ax0
ax = axs[0]
taxonomy_draw(ax, refined, "fPos_y", "fNeg_y", **params)
# for i in np.arange(0, 0.25, 0.025):
#     ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)
# ax1
ax = axs[1]
taxonomy_draw(ax, refined, "fPos_x", "fNeg_x", **params)
# for i in np.arange(0, 0.25, 0.025):
#     ax.plot([0, 0.3 - i], [i ,0.3], linestyle = "--", color= "k", zorder = 5)


# %% codecell
top_gp = refined.phylum.value_counts().index[:4]
refined[refined.phylum.isin(top_gp)].groupby("phylum").mean()[["NCD1000", "ncd1000", "MWkDa", "mwkda"]]
