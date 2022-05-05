#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import HALO_DIR, TARGETS
from src.search.data import load_proteomedata, taxcols
from src.search import fasta

# %% codecell
proteomes = load_proteomedata("uniprot").drop(taxcols, axis = 1)
hits = pd.read_csv(HALO_DIR / "bacteria/uniprot/iter1.csv")
hits["mwkda"] = hits.mass / 1000
hits["ncd1000"] = hits.ncd * 1000
hits = hits.sort_values("evalue") \
           .drop_duplicates(["proteomecode"], keep = "first") \
           .reset_index(drop=True)


# %% codecell
hits.describe()


# %% codecell
pattern = "cyanin|azurin|cupredoxin|cyanin"
hits["choice"] = hits.name.str.contains(pattern, case = False)


# %% codecell
hits[(hits.choice)].mwkda.plot.kde()


# %% codecell
hits[(hits.choice)].evalue.plot.kde().set_xscale("log")


# %% codecell
sns.scatterplot(data=hits[hits.choice], x="ncd1000", y="mwkda")


# %% codecell
fasta.write_hits(hits[hits.choice].to_dict("records"), TARGETS / "halocyanin_bacteria_uniprot.fasta")
