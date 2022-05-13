#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import HALO_DIR, PROSITES, ROOT
from src.search import prosite, filter, fasta

# %% codecell
# hits in iter 1
with open(HALO_DIR / "archaea/uniprot/iter1.csv", "r") as csvfile:
    hits = [row for row in csv.DictReader(csvfile)]
    hits = pd.DataFrame(hits)

# models
TEMP = ROOT / "analysis/searches/halocyanin/temp/"
p = TEMP / "halocyanin_models_idlist"
modelids = [
        line.strip()
        for line in open(p, "r").readlines()
    ]
# prosites
ps00078 = prosite.parse(PROSITES / "ps00078_halocyanin_archaea_uniprot.csv")
ps00196 = prosite.parse(PROSITES / "ps00196_halocyanin_archaea_uniprot.csv")


# %% codecell
# Hits, models and motifs
modelled = len(modelids)
checked = np.intersect1d(modelids, hits.id)
models_w_cua = np.intersect1d(checked, ps00078.protein)
models_w_type1 = np.intersect1d(checked, ps00196.protein)
hits_w_type1 = np.intersect1d(hits.id, ps00196.protein)
hits_w_type1_n_cua = np.intersect1d(hits_w_type1, ps00078.protein)


print(
    f"""
    Gabrielle made AF2 models for {modelled} proteins.
    {checked.shape[0]} were checked for prosite motifs:
        {models_w_cua.shape[0]} contained a CuA center,
        {models_w_type1.shape[0]} contained a Type1 copper center,

    Had we selected for the presence of the Type1 motif in our targets
    we would have obtained:
        a sample of {hits_w_type1.shape[0]} sequences.
        only {hits_w_type1_n_cua.shape[0]} of those sequences also contain the CuA motif.
    """
)


# %%codecell
wouldbehits_w_type1 = hits.set_index("id") \
                          .reindex(ps00196.protein)\
                          .sort_values("evalue")\
                          .drop_duplicates("proteomecode", keep="first")\
                          .index

wouldbehits_w_type1_n_cua = np.intersect1d(wouldbehits_w_type1, ps00078.protein)


# %% codecell
# Write fasta with hits contain both motifs for AF2 modelling
cua_type1 = hits[hits.id.isin(hits_w_type1_n_cua)]
# fasta.write_hits(cua_type1.to_dict("records"), TEMP / "cua_type1.fasta")


# %% codecell
ps = pd.concat([ps00078, ps00196])
hits.merge(ps, left_on="id", right_on="protein", how="inner")\
    .set_index("id").loc[hits_w_type1_n_cua][["motif", "hit_start", "hit end"]]



############################# CLUSTER II #######################################
# %% codecell
cii = pd.read_csv(HALO_DIR / "clusterii/ncbi/iter1.csv" )
ps00196 = prosite.parse(PROSITES / "ps00196_halocyanin_clusterii_ncbi.csv")

cii.merge(ps00196, left_on="id", right_on="protein", how="inner")\
   .sort_values("evalue")\
   .drop_duplicates("proteomecode", keep = "first").name.value_counts()
    
