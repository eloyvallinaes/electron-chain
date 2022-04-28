"""
HMMER searching with pyhmmer
"""

import pyhmmer
import pandas as pd
from src.config import DATA_DIR, QUERY_DIR, ROOT

alphabet = pyhmmer.easel.Alphabet.amino()
background = pyhmmer.plan7.Background(alphabet)
builder = pyhmmer.plan7.Builder(alphabet)


def load_sequence(filepath: str):
    """
    Load sequences from file into DigitalSequence formats.

    :param filepath: a path to fasta file containing as least one sequence
    :return: a list of DigitalSequence objects (pyhmmer)
    """
    with pyhmmer.easel.SequenceFile(filepath, alphabet=alphabet) as sf:
        sequences = list(sf)

    return [s.digitize(alphabet) for s in sequences]


def load_msa(filepath: str):
    """
    Load MSA from file in stockholm format.

    :param filepath: a path to stockholm file containing the alignment
    :return: pyhmmer.easel.MSA object
    """
    with pyhmmer.easel.MSAFile( alignfile,
                                format = "stockholm",
                                digital=True,
                                alphabet=alphabet) as msa_file:
        msa = next(msa_file)
        msa.name = b"placeholder" # must be set

    return msa


def parse_hit(hit, db):
    if db == "uniprot":
        return dict(id=str(hit.name).split("|")[1], evalue=hit.evalue)
    elif db == "ncbi":
        return dict(id=str(hit.name).split("_")[2], evalue=hit.evalue)


def runsearch(mode: str, query: list, targets: list, db: str):
    """
    Wrapper for pyhmmer.hmmer.phmmer and pyhmmer.hmmer.hmmsearch depending on
    mode parameter.

    :param mode: search strategy: 'protein' for phmmer; 'msa' for hmmsearch
    :param query: a list of len 1 containing:
        a DigitalSequence object (mode = 'protein');
        a DigitalMSA object (mode = 'msa')
    :param targets: a list of len > 1 containing DigitalSequences
    :return: a list of dictionaries with id and evalue fields
    """
    if len(query) > 1 or len(targets) <= 1:
        message = "Too many queries or too few targets"
        raise ValueError(message)

    if mode == "protein":
        hits = list( pyhmmer.hmmer.phmmer(query, targets) )[0]
    elif mode == "msa":
        hits = list( pyhmmer.hmmer.hmmsearch(query, targets) )[0]
    else:
        ValueError(f"mode must be either 'protein' or 'msa' not {mode}")

    results = []
    for hit in hits:
        results.append(parse_hit(hit, db))
    return results


def main(mode, queryfilepath, db, group, outdir, n=0):
    """
    Top searching interface for applying a search strategy (mode) with a
    query (queryfilepath). Expects targets in directory DATA_DIR / db / group
    and writes outputs to directory outdir. For testing shorter searches use
    keyword n.
    """
    if mode == "protein":
        query = load_sequence(queryfilepath)
    elif mode == "msa":
        query = load_msa(queryfilepath)
    else:
        ValueError(f"mode must be either 'protein' or 'msa' not {mode}")

    targets_dir = DATA_DIR / db / group
    results = []
    c = 0
    for file in targets_dir.iterdir():
        print(f"Searching in: {file.stem}")
        targets = load_sequence(file)
        hits = runsearch(mode, query, targets, db)
        for hit in hits:
            hit.update( {"proteomecode" : file.stem} )

        results.extend(hits)
        c +=1
        if n > 0 and c == n:
            break

    return results



if __name__ == '__main__':
    # test run
    mode = "protein"
    db = "uniprot"
    group = "Archaea"
    query = QUERY_DIR / "P39442.fasta"
    outdir = ROOT / 'searches/halocyanin/archaea/uniprot'
    n=10
    results = main(mode, query, db, group, outdir, n = n)
    # write
    pd.DataFrame(results).to_csv(outdir / f"test{n}.csv", index=False)
