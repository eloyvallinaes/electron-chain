#!/usr/bin/python
"""
Apply selection rules to list of hits.

"""

import numpy as np
import pandas as pd

def top_per_proteome(hits):
    """
    Sort hist by evalue and keep lowest scoring hit per proteome.

    :param hits: list of dictionaries containing at least the fields evalue and
                 proteomecode
    :returns: pd.DataFrame

    """
    hitsdf = pd.DataFrame(hits)
    return hitsdf.sort_values("evalue", ascending=True) \
                 .drop_duplicates(["proteomecode"], keep="first")


def main(hits, pattern=".*", mw=np.Inf, n=0):
    """
    Sort hist by evalue and keep lowest scoring hit per proteome. Optionally
    match protein names to a pattern, apply a molecular weight threshold or
    limit result to n entries.

    :param hits: list of dictionaries containing at least evalue, proteomecode,
                  name and mw columns.
    :param pattern: a string pattern to match protein names to.
    :param mw: a numeric value to use as cap for protein mw column.
    :param n: retrieve at most n rows.
    :returns: list of dictionaries
    """
    top = top_per_proteome(hits)
    selector = ( top.name.str.contains(pattern) ) & (top.mass < mw)
    top = top[selector]
    if n >0:
        top = top.head(n)
    return top.to_dict('records')
