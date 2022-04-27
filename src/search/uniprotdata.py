#!/usr/bin/python
"""
Query UNIPROT API with protein accessions and return CSV file with
results.

Usage: python3 uniprotdata infilepath outfilepath

Arguments:
    infilepath - csv file with parsed results from HMMER searches
    outfilepath - csv file with uniprotdata and sequences

"""
import os
import sys
import csv
import requests
import pandas as pd

KBCOLUMNS = [
                'id',
                'entry name',
                'mass',
                'sequence',
                'protein names',
                'lineage(SUPERKINGDOM)',
                'lineage(PHYLUM)',
                'lineage(CLASS)',
                'lineage(ORDER)',
                'lineage(FAMILY)',
                'lineage(GENUS)',
                'lineage(SPECIES)',
        ]

KEYS = [
        'upacc',
        'entry',
        'mass',
        'sequence',
        'name',
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'id'
    ]

def parse(r):
    reader = csv.reader(  r.text.splitlines(), delimiter = "\t" )
    next(reader)
    records = [
                {
                    key: val
                    for key,val in zip(KEYS, line)
                }
                for line in reader
            ]
    for record in records:
        record['mass'] = float( record['mass'].replace(",", "") )
    return records


def retrieve(accs: list, db):
    url = 'https://www.uniprot.org/uploadlists/'

    if db not in ["uniprot", "ncbi"]:
        raise ValueError(f"db must be either 'ncbi' or 'uniprot' not {db}")

    idfrom = "ACC" if db == "uniprot" else "EMBL"
    url = 'https://www.uniprot.org/uploadlists/'
    query = " ".join(accs)
    params = {
                'from': idfrom, # "ACC" for uniprot, "EMBL" for ncbi
                'to': 'ACC',
                'format': 'tab',
                'query': query,
                'columns':','.join(KBCOLUMNS),
    }
    r = requests.post(url, data=params)
    records = parse(r)
    return records


def main(infilepath, outfilepath, db):
    if dbin not in ["uniprot", "ncbi"]:
        raise ValueError(f"db must be either 'ncbi' or 'uniprot' not {db}")

    dtype = {"mass": float}
    idfrom = "ACC" if db == "uniprot" else "EMBL"
    with open(infilepath) as csvfile:
        accs = [line.split(",")[0] for line in csvfile]
        records = retrieve(accs, idfrom)
        pd.DataFrame(records).astype(dtype) \
          .to_csv(outfilepath, index = False)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3]) # infilepath, outfilepath, idfrom
