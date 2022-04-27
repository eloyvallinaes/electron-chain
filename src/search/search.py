"""
HMMER searching with pyhmmer
"""

import pyhmmer
import pandas as pd
from src.config import DATA_DIR, QUERY_DIR, ROOT

alphabet = pyhmmer.easel.Alphabet.amino()


def load_sequence(filepath: str):
    """
    Load sequences from file into DigitalSequence formats.

    :param filepath: a path to fasta file containing as least one sequence
    :return: a list of DigitalSequence objects (pyhmmer)
    """
    with pyhmmer.easel.SequenceFile(filepath, alphabet=alphabet) as sf:
        sequences = list(sf)

    return [s.digitize(alphabet) for s in sequences]


def psearch(query: list, targets: list, db: str):
    """
    Wrapper for pymmer.hmmer.phmmer.
    :param query: a list of len 1 containing a DigitalSequence objects
    :param targets: a list of len > 1 containing DigitalSequences
    :return: a list of dictionaries with id and evalue fields
    """
    if len(query) > 1 or len(targets) <= 1:
        message = "Too many queries or too few targets"
        raise NotImplementedError(message)

    hits = list( pyhmmer.hmmer.phmmer(query, targets) )[0]

    results = []
    for hit in hits:
        results.append(parse_hit(hit, db))
    return results

def parse_hit(hit, db):
    if db == "uniprot":
        return dict(id=str(hit.name).split("|")[1], evalue=hit.evalue)
    elif db == "ncbi":
        return dict(id=str(hit.name).split("_")[2], evalue=hit.evalue)


def main(queryfilepath, db, group, outdir, n=0):
    query = load_sequence(queryfilepath)
    targets_dir = DATA_DIR / db / group
    results = []
    c = 0
    for file in targets_dir.iterdir():
        if c==n and n >0:
            break
            
        print(file.stem)
        targets = load_sequence(file)
        hits = psearch(query, targets, db)
        for hit in hits:
            hit.update( {"proteomecode" : file.stem} )

        results.extend(hits)
        c +=1

    return results


if __name__ == '__main__':
    db = "uniprot"
    group = "Archaea"
    query = QUERY_DIR / "P39442.fasta"
    outdir = ROOT / 'searches/halocyanin/archaea/uniprot'
    results = main(query, db, group, outdir)
    write_results(results, outdir / "iter0.csv")
