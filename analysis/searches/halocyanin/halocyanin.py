#!/usr/bin/python

import os
import csv
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
DIR = "/Users/eloyvallina/Documents/Dropbox/phd/Projects/genomescreening"
from analysis.data import load_carrier_data

# %% codecell
up = pd.read_csv("csv/halocyanin/Archaea/iter1_up.csv").set_index("acc")
info = pd.read_csv("csv/halocyanin/Archaea/iter1.csv").set_index("acc")

df = pd.concat([up, info], axis = "columns")

# %% codecell
df.head()
dff = df.sort_values("evalue").drop_duplicates(subset=["upcode"], keep = "first")
dff[dff["order"] == "Candidatus Poseidoniia"]


# %% codecell
top = [ "Halobacteria", "Nanohaloarchaea", "Candidatus Poseidoniia",
        "Methanomicrobia", "Thermoplasmata", 'nan']
cc = [  "tab:blue", "fuchsia", "tab:red", "tab:green", "tab:orange", "black"]
colors = {name: c for name, c in zip(top, cc)}

def draw(df, x, y, category):
    fig = plt.figure(figsize=(6,6))
    ax  = fig.subplots(1,1)
    for name, subset in df.groupby(category, dropna=False):
        name = str(name)
        if name in top:
            ax.scatter( subset[x], subset[y], label = name, color = colors[name],
                        edgecolor = "w", marker = "o", linewidth = 0.2)
        else:
            ax.scatter( subset[x], subset[y], label = name, color = "gray",
                        edgecolor = "w", marker = "o", linewidth = 0.2)
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_xlabel(x)
    ax.set_ylabel(y)

    return fig

# %% codecell
selector = ( df.name.str.contains("[Cc]yanin") ) | \
           ( df.name.str.contains("[Uu]ncharacterized") ) | \
           ( df.name.str.contains("[Cc]opper") ) | \
           ( df.name.str.contains("[Cc]upredoxin") ) | \
           ( df.name.str.contains("[Aa]zurin") )

filter =   ( df.name.str.contains("[Tt]ransport") ) | \
           ( df.name.str.contains("[Oo]xidase") ) | \
           ( df.name.str.contains("[Kk]azal") )

# %%codecell
fig = draw(df[~filter], x = "ncd1000", y="mwkda", category="class")
fig.axes[0].set_title("Protein values")

# %%codecell
fig = draw(df[~filter], x = "fCharged_x", y="fFatty_x", category="class")
fig.axes[0].set_title("Protein values")

#%% codecell
fig = draw(df[~filter], x = "NCD1000", y="MWkDa", category="class")
fig.axes[0].set_title("Avg. proteome values")

# %% codecell
df[~filter].acc.values

# %% codecell
def writeFasta(accs, outname, upname):
    with open(outname, "w") as f:
        with open(upname, "r") as sequences:
            reader = csv.DictReader(sequences)
            for line in reader:
                if line["acc"] in accs:
                    acc = line["acc"]
                    entry = line["entry"]
                    mass = float(line["mass"]) / 1000
                    name = line["name"]
                    seq = line["sequence"]
                    f.write(f">{acc}|{entry}|MW={mass}kDa|{name}\n{seq}\n")

# %% codecell
# writeFasta(df[~filter].acc.values, "fasta/halocyanin_filtered.fasta", "csv/iter2_up.csv")

df.shape
