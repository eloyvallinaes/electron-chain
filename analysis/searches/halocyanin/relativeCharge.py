#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from src.config import HALO_DIR, PROSITES
from src.search import prosite


# %% codecell
hits = pd.read_csv(HALO_DIR / "archaea/uniprot/iter1.csv" )
prosites = prosite.parse(PROSITES / "ps00078_halocyanin_archaea_uniprot.csv")
hits = hits.merge(prosites, left_on="id", right_on="protein", how="left")

prosites.duplicated(subset="protein").sum()
