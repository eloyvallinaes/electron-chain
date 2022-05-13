"""
Characterisation of targets selected by Nannan for crystalisation and
experimental research.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.config import ROOT
from src.search.data import load_proteomedata, isocharge_draw, taxcols

# paths
HALO = ROOT / "searches" / "halocyanin"
TEMP = ROOT / "analysis" / "searches" / "halocyanin" / "temp"
# proteomes
proteomes = load_proteomedata("uniprot").drop(taxcols, axis = 1)
# targets
up = pd.read_csv(TEMP / "potential_c1a_halocyanin_targets_up.csv")
phys = pd.read_csv(TEMP / "potential_c1a_halocyanin_targets_phys.csv")
phys["ncd1000"] = phys.ncd * 1000
phys["mwkda"] = phys.mass / 1000
targets = up.merge(phys, on = "upacc")
# combine
targets = targets.merge(proteomes, left_on="upproteome", right_on="Name", how = "left")

# %% codecell
fig = plt.figure(figsize=(6,6))
ax = fig.subplots(1,1)
isocharge_draw(ax, targets)


# %% codecell
cols =  [   "upacc", "ncd1000", "NCD1000", "mwkda", "MWkDa", "GC", "fNeg_x",
            "fNeg_y", "fPos_x", "fPos_y", "species"
    ]
names = [   "upacc", "protein_ncd1000", "proteome_ncd1000", "protein_mwkda",
            "proteome_mwkda", "GC", "fNeg_protein", "fNeg_proteome",
            "fPos_protein", "fPos_proteome", "species"
    ]
table = targets[cols]
table.columns = names
table.to_csv("potential_c1a_halocyanin_targets_properties.csv", index=False)
