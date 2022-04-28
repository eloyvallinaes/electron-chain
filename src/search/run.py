"""
Top level script to run a search pipeline
"""

import pandas as pd
from src.config import QUERY_DIR, ROOT
from src.search import search, uniprotdata, measure, filter, fasta, clustalo

def hit_uniprot_update(hits):
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
    db = "uniprot"
    group = "Archaea"
    query = QUERY_DIR / "P39442.fasta"
    outdir = ROOT / 'searches/halocyanin/archaea/uniprot'

    # ITER 0
    # iter0 = search.main("protein", query, db, group, outdir)
    # hit_uniprot_update(iter0)
    # hit_measure_update(iter0)
    # write_outputs(iter0, outdir, "iter0", alignment=True)

    # ITER 1
    msafile = outdir / "iter0.sto"
    iter1 = search.main("msa", msafile, db, group, outdir)
    hit_uniprot_update(iter1)
    hit_measure_update(iter1)
    write_outputs(iter1, outdir, "iter1", alignment=False)
