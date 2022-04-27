#!/usr/bin/python
"""
Apply selection rules to list of hits.

"""

import pandas as pd

def top_per_proteome(hits):
    hitsdf = pd.DataFrame(hits)
    return hitsdf.sort_values("evalue", ascending=True) \
                 .drop_duplicates(["proteomecode"], keep="first")


def main(hits, pattern, mw):
    top = top_per_proteome(hits)
    selector = ( top.name.str.contains(pattern) ) & (top.mass < mw)
    return top[selector].to_dict('records')
