#!/usr/bin/python
import csv
import sys
import pandas as pd

class Protein():
    def __init__(self, upacc, id, entry, mass, sequence, name):
        self.upacc = upacc
        self.id = id
        self.entry = entry
        self.name = name
        self.mass = float(mass)
        self.sequence = sequence
        self.length = len(self.sequence)
        self.sasa = self._sasa()
        self.pos = self.sequence.count("K") + self.sequence.count("R")
        self.neg = self.sequence.count("D") + self.sequence.count("E")
        self.net = self.pos - self.neg
        self.ncd = self.net / self.sasa
        self.fat = sum([self.sequence.count(letter) for letter in "FLIV"])

    def _sasa(self):
        """Miller's approx"""
        return 6.3 * self.mass**0.73

    def serialize(self):
        return {
                "upacc": self.upacc,
                "id" : self.id,
                "entry": self.entry,
                "name": self.name,
                "mass": self.mass,
                "sasa": self.sasa,
                "fPos": self.pos / self.length,
                "fNeg": self.neg / self.length,
                "net" : self.net,
                "ncd" : self.ncd,
                "fFatty": self.fat / self.length,
            }


def measurement(record: dict):
    keys = ['upacc', 'id', 'entry', 'mass', 'sequence', 'name']
    r = {key: record[key] for key in keys}
    p = Protein(**r)
    return p.serialize()


def measure_file(csvfile):
    keys = ['upacc', 'id', 'entry', 'mass', 'sequence', 'name']
    phys = []
    reader = csv.DictReader(open(csvfile, "r"))
    for record in reader:
        r = {key: record[key] for key in keys}
        p = Protein(**r)
        data = p.serialize()
        phys.append(data)
    return phys


def main(infilepath, outfilepath):
    phys = measure(infilepath)
    pd.DataFrame(phys).to_csv(outfilepath, index=False)


if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])
