#!/usr/bin/python
"""
An interface for obtaining basic clustal omega alignments through
the EBI REST API.
"""

import pyhmmer
import requests
from time import sleep
from src.config import UP_DIR, HALO_DIR

URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/"
EMAIL = "eloyvallina33@gmail.com"


def submit(fastafile):
    data = {
            "sequence": open(fastafile, "r").read(),
            "email": EMAIL,
            "outfmt": "stockholm",
        }
    r = requests.post(URL + "run/", data=data)
    if not r.ok:
        print(data)
        raise Exception(f"Request failed: {r.text}")

    jobid = r.content.decode()
    return jobid


def check_finished(jobid):
    r = requests.get(URL + f"status/{jobid}")
    if not r.ok:
        raise Exception(f"Request failed: {r.text}")

    elif r.content.decode() == "FINISHED":
        return True

    return False


def retrieve(jobid, stofilepath):
    r = requests.get(URL + f"result/{jobid}/aln-stockholm")
    if not r.ok:
        raise Exception(f"Request failed: {r.text}")

    content = r.content.decode()
    with open(stofilepath, "w") as stofile:
        stofile.write(content)


def main(fastafile, stofilepath):
    jobid = submit(fastafile)
    while not check_finished(jobid):
        sleep(30)

    retrieve(jobid, stofilepath)




if __name__ == '__main__':
    fastafile = HALO_DIR / "archaea/uniprot/iter0.fasta"
    alignfile = HALO_DIR / "archaea/uniprot/iter0.sto"
    main(fastafile, alignfile)
