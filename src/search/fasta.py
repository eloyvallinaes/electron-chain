#!/usr/bin/python

"""
Functions related to sequences in FASTA format
"""
from Bio import SeqIO


def clean_character(filepath, character):
    """
    Remove sequence records containing a character and overrite file.
    Intended to clean up "-" characters in CDS from NCBI Assemblies
    """
    records = SeqIO.parse(filepath, "fasta")
    passed = [r for r in records if not character in r.seq]
    SeqIO.write(passed, filepath, "fasta")


def write_hits(hits: list, filepath):
    """
    Format hits list-of-dicts to FASTA and write file. hits must contain fields:
    id, name, mass and sequence
    """
    with open(filepath, "w") as sequences:
        for hit in hits:
            acc = hit['id']
            name = hit['name']
            mw = hit['mass']
            seq = hit['sequence']
            sequences.write(f">{acc}|{name}|MW={mw}Da\n{seq}\n")
