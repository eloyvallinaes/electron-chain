"""
Top level script to run a search pipeline
"""

import pandas as pd
from src.config import QUERY_DIR, ROOT
from src.search import search, uniprotdata, measure, filter, fasta

db = "ncbi"
group = "ClusterII"
query = QUERY_DIR / "P39442.fasta"
outdir = ROOT / 'searches/halocyanin/clusterii/ncbi'

# first search
iter0 = search.main(query, db, group, outdir)

# map hits to uniprot
accs = [hit['id'] for hit in iter0]
records = uniprotdata.retrieve(accs, db)
for hit in iter0.copy():
    for record in records:
        if hit['id'] == record['id']:
            hit.update(record)
            break
    else:
        iter0.remove(hit) # no up mapping found

# measure physicochemical properties
for hit in iter0:
    try:
        phys = measure.measurement(hit)
        hit.update(phys)
    except KeyError:
        print(f"Couldn't measure {hit}")


# selection rules: top hit per proteome, name contains pattern and mw threshold
# iter0 = filter.main(iter0, "[Hh]alocyanin", 50000)

# write output
pd.DataFrame(iter0).to_csv(outdir / "iter0.csv", index=False) # csv file
fasta.write_hits(iter0, outdir / "iter0.fasta")
