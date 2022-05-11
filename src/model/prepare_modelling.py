import os
import sys
import argparse
from src.config_gp import *
from src.model.manager import FastaManager

def main():
    p = argparse.ArgumentParser(
        description = '- Analyze fasta and setup folder for modelling -')
    p.add_argument('-i', required=True, help='input multi-sequence fasta')
    ns = p.parse_args()

    manager = FastaManager(ns.i, ROOT+'/models/')
    manager.run_clustering(MMSEQ)
    manager.pick_fasta()

if __name__ == '__main__':
    main()

