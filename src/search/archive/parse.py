#!/usr/bin/python
"""
Parse hmmer result TSV files to a single CSV file listing
upcode, taxid, accession and evalue.

Usage:
python3 parse.py directory outname

Arguments:
    directory - the location of the tsv files produced by hmmer
    outfilepath - parsed results with .csv extension
"""
import os
import sys
import pandas as pd
from Bio import SearchIO

def find_files(path):
    return os.listdir(path)


def parse_hit(hit, dbin):
    if dbin == "uniprot":
        return dict(id=hit.id.split("|")[1], evalue=hit.evalue)
    elif dbin == "ncbi":
        return dict(id=hit.id.split("_")[2], evalue=hit.evalue)


def parse_filename(filename, dbin):
    filename = filename.replace(".tsv", "")
    if dbin == "uniprot":
        upcode, taxid = filename.split("_")
        return dict(proteomecode=upcode, taxid=taxid)
    elif dbin == "ncbi":
        return dict(proteomecode=filename, taxid="")


def parse_searches(directory, dbin):
    records = []
    for name in find_files(directory):
        data = parse_filename(name, dbin)
        path = os.path.join(directory, name)
        for result in SearchIO.parse(path, format = "hmmer3-tab"):
            for hit in result.hits:
                record = parse_hit(hit, dbin)
                record.update(data)
                records.append(record)
    return records


def main(directory, outfilepath, dbin):
    if dbin not in ["ncbi", "uniprot"]:
        raise ValueError(f"dbin must be either 'ncbi' or 'uniprot' not {dbin}")

    records = parse_searches(directory, dbin)
    pd.DataFrame(records).to_csv(outfilepath, index = False)

if __name__ == '__main__':
    main( sys.argv[1], sys.argv[2], sys.argv[3] ) # directory, outfilepath, dbin
