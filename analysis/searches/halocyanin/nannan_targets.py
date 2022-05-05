"""
Characterisation of targets selected by Nannan for crystalisation and
experimental research.
"""

import numpy as np
import pandas as pd
from src.config import ROOT
from src.search.data import load_proteomedata, isocharge_draw

HALO = ROOT / "searches" / "halocyanin"
TARGETS = ROOT / "analysis" / "searches" / "halocyanin" / "halocyanin_targets.csv"


# %% codecell
up = load_proteomedata("uniprot")
ncbi = load_proteomedata("ncbi")
ncbi["Name"] = ncbi.AssemblyID
proteomes = pd.concat([up, ncbi])


# %% codecell
a0 = pd.read_csv(HALO / "archaea" / "uniprot" / "iter1.csv", header = 0)
a1 = pd.read_csv(HALO / "clusterii" / "ncbi" / "iter0.csv", header = 0)
h = pd.concat([a0, a1])
nannan = pd.read_csv(TARGETS, header = 0)


#%% codecell
taxcols = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
targets = nannan.merge(h, on="upacc", how = "left") \
                .merge(
                        proteomes.drop(taxcols, axis = 1),
                        left_on="proteomecode",
                        right_on="Name",
                        how = "left",
                    )
targets[targets.group == "cluster1b"][["species", "upacc"]]


# %% codecell
f = isocharge_draw(targets, "kingdom", "Archaea", "group")
f.axes[0].legend(bbox_to_anchor=(1,1))
