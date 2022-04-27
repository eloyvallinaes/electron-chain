#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
DIR = "/Users/eloyvallina/Documents/Dropbox/phd/Projects/genomescreening"
sys.path.append(os.path.join(DIR, "utils"))
from data import load_uniprot_data, load_ncbi_data

def load_proteins(target: str, group: str, i: int):
    dtype={'taxid': str}
    dir = os.path.join("csv", target, group)
    infopath = os.path.join(dir, f"iter{str(i)}.csv")
    physpath = os.path.join(dir, f"iter{str(i)}_phys.csv")
    info = pd.read_csv(infopath, header=0, dtype=dtype)
    phys = pd.read_csv(physpath, header=0)
    # extra columns
    phys["ncd1000"] = phys.ncd * 1000
    phys["mwkda"] = phys.mass / 1000
    phys["fCharged"] = phys.fPos + phys.fNeg
    return phys.merge(info, on = "id", how = "left")

def load_carrier_data(target: str, group: str, i: int, db: str):
    """
    Produce merged dataframe of physicochemical coordinates from proteins
    and proteomes.

    Each UP code is represented by its top hit (evalue)
    """
    # Proteome data
    if db not in ["uniprot", "ncbi"]:
        raise ValueError(f"db must be 'uniprot' or 'ncbi' not {db}")

    proteomes = load_uniprot_data() if db == "uniprot" else load_ncbi_data()

    # Protein data
    first = load_proteins(target, group, i)

    # Merge
    second = first.merge(proteomes, left_on = "proteomecode", how="inner", right_on="Name")
    # reduce to top hit per UP code by evalue
    second.sort_values("evalue", inplace=True)
    second.drop_duplicates("proteomecode", keep = "first", inplace=True)
    return second.reset_index(drop=True)
