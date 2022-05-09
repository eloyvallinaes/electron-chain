"""
Top level script to run a search pipeline
"""

import os
import json
import pathlib
import argparse
import numpy as np
import pandas as pd
from src.config import QUERY_DIR, ROOT
from src.search import search, uniprotdata, measure, filter, fasta, clustalo

def hit_uniprot_update(hits, db):
    """
    Map proteins in hits to uniprot and update list records
    """
    accs = [hit['id'] for hit in hits]
    records = uniprotdata.retrieve(accs, db)
    for hit in hits.copy():
        for record in records:
            if hit['id'] == record['id']:
                hit.update(record)
                break
        else:
            hits.remove(hit) # no up mapping found

def hit_measure_update(hits):
    """
    Measure physicochemical properties of sequences and update list records
    """
    for hit in hits:
        phys = measure.measurement(hit)
        hit.update(phys)

def write_outputs(hits, outdir, name, alignment=False):
    pd.DataFrame(hits).to_csv(outdir / f"{name}.csv", index=False) # csv file
    fasta.write_hits(hits, outdir / f"{name}.fasta")
    if alignment:
        print("Aligning hits")
        clustalo.main(outdir / f"{name}.fasta", outdir / f"{name}.sto")


if __name__ == '__main__':
    description = "Run a search pipeline mimicking jackhmmer."
    parser = argparse.ArgumentParser(
        description="Run a search pipeline"
    )
    parser.add_argument(
        'db',
        choices=["uniprot", "ncbi"]
    )
    parser.add_argument(
        'group',
        type=str,
        help="A subset of db where to search"
    )
    parser.add_argument(
        'query',
        type=pathlib.Path,
        help="Fasta file with single sequence to initialise the search"
    )
    parser.add_argument(
        'outdir',
        type=pathlib.Path,
        help="A directory where to store outputs"
    )
    parser.add_argument(
        '--resume',
        action='store_true',
        help="Resume iter 1 assuming iter 0 alignment is present in outdir",
    )
    parser.add_argument(
        '--filter-name',
        default = ".*",
        help = "A pattern for name filtering after iter 0., Default is '.*'"
    )
    parser.add_argument(
        '--filter-mw',
        type=float,
        default=np.inf,
        help="A MW threshold in Da to filter hits after iter 0. Default is np.inf"
    )

    args = parser.parse_args()

    # CREATE OUTPUT DIR IF NECESSARY
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # WRITE LOG
    with open(args.outdir / "log.txt", "w") as logfile:
        logfile.write("Ran search.run with parameters:")
        for k, v in vars(args).items():
            logfile.write(f'{k}: {v}\n')

    if not args.resume:
        # ITER 0
        iter0 = search.main("protein", args.query, args.db, args.group, args.outdir)
        hit_uniprot_update(iter0, db=args.db)
        hit_measure_update(iter0)
        iter0 = filter.main(
                            iter0,
                            args.filter_name,
                            args.filter_mw,
                            4000 # max sequences allowed by clustalw
                        )
        write_outputs(iter0, args.outdir, "iter0", alignment=True)
        print("ITER 0: Done")

    # ITER 1
    msafile = args.outdir / "iter0.sto"
    iter1 = search.main("msa", msafile, args.db, args.group, args.outdir)
    hit_uniprot_update(iter1, db=args.db)
    hit_measure_update(iter1)
    write_outputs(iter1, args.outdir, "iter1", alignment=False)
    print("ITER 1: Done")
