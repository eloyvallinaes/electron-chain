#!/usr/bin/python
"""
Basic parsing of prosite results in table format.
See https://prosite.expasy.org/scanprosite/)
"""

import pandas as pd
from src.config import PROSITES

columns = [
        'protein',
        'hit_start',
        'hit end',
        'motif',
        'scorex',
        'scorey',
        'confidence',
        'region'
    ]

def parse(filepath):
    df = pd.read_csv(filepath, delimiter='\s+', names=columns)
    df['protein'] = df.protein.str.split('|', expand=True)[0]
    df['confidence'] = df.confidence.map({"(0)": 'high', "(-1)": 'low'})
    return df


if __name__ == "__main__":
    print(parse(PROSITES / "ps00078_halocyanin_archaea_uniprot.csv"))
