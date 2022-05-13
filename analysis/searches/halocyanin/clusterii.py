#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.search import prosite
from src.config import HALO_DIR, PROSITES
from src.search.data import load_proteomedata, taxcols, isocharge_draw, taxonomy_draw


# %% codecell
targets = pd.read_csv(HALO_DIR / "clusterii/ncbi/iter1.csv")
targets["mwkda"] = targets.mass / 1000
targets["ncd1000"] = targets.ncd * 1000
ncbi = load_proteomedata("ncbi")
ps00196 = prosite.parse(PROSITES / "ps00196_halocyanin_clusterii_ncbi.csv")
ps00196.head()

# %% codecell
def target_proteome_merge(targets, proteomes):
    return targets.merge(
                    proteomes.drop(taxcols, axis = 1),
                    left_on="proteomecode",
                    right_on="Name",
                    how = "left"
                )


# %% codecell
# List for Nannan
columns = [
        'id', 'ncd1000', 'mwkda', "MWkDa", "NCD1000", "class", "species",
        "upacc", "name"
    ]
selection = ( df["class"] == "Candidatus Poseidoniia" ) & ( df.mwkda < 40 )
narrow = df[selection]
narrow[columns]
# .to_csv("cii_halocyanin_targets.csv", index=False)


# %% codecell
# 5 targets selected for experiments by Nannan
accs = ["A0A6V8ERJ7", "A0A822XG38", "A0A7N4A687", "A0A6V8ELI7", "A0A3A5V768"]
narrow = narrow.merge(ps00196, left_on="id", right_on="protein", how = "left")
narrow.set_index("upacc").reindex(accs)

# %%codecell
narrow.motif = narrow.motif.fillna("None")


# %% codecell
fig = plt.figure()
ax = fig.subplots(1,1)
isocharge_draw(ax, narrow, lower = "motif")
ax.legend()
