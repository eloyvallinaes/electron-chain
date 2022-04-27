"""
Characterisation of targets selected by Nannan for crystalisation and
experimental research.
"""

import pandas as pd
from src.config import ROOT
from src.search import data

HALO = ROOT / "searches" / "halocyanin"
TARGETS = ROOT / "analysis" / "searches" / "halocyanin" / "halocyanin_targets.csv"

# %% codecell
h0 = pd.read_csv(HALO / "archaea" / "uniprot" / "iter0.csv", header = 0)
h1 = pd.read_csv(HALO / "clusterii" / "ncbi" / "iter0.csv", header = 0)
targets = pd.read_csv(TARGETS, header = 0)


#%% codecell
targets.head()
h0.shape
h1.shape
