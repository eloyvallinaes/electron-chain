#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import HALO_DIR, CYTO_DIR
from src.search.data import load_proteomedata, isocharge_draw, taxcols

# %% codecell
cytob = pd.read_csv(CYTO_DIR / "bacteria" / "uniprot" / "iter1.csv") \
          .sort_values("evalue", ascending=True) \
          .drop_duplicates(["proteomecode"], keep="first") \
          .reset_index()


# %% codecell
up = load_proteomedata("uniprot").drop(taxcols, axis = "columns")
df = cytob.merge(up, left_on="proteomecode", right_on="Name", how="inner")


# %% codecell
pattern = "(?<!bc1 complex )cytochrome c (?!oxidase|reductase|oxidoreductase|peroxidase)"
phyla = cytob.phylum.value_counts()[:10].index
selector = ( cytob.name.str.contains(pattern, case=False) ) & \
           ( cytob.evalue < 0.001 ) & \
           ( cytob.phylum.isin(phyla) )

subset = cytob[selector].name.unique()


# %% codecell
fig = plt.figure(figsize=(6,6))
ax  = fig.subplots(1,1)
sns.scatterplot(data=subset, x="ncd", y="mass", hue="phylum")

# %% codecell
fig = isocharge_draw(df.sample(100), "kingdom", "Bacteria", "phylum", n=20)
