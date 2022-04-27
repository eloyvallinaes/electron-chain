#!/usr/bin/python
"""
Download ref_proteomes from Uniprot FTP
"""

import requests
from src.config import UP_DIR

URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"

def get_list(kingdom):
    url = f"{URL}/README"
    r = requests.get(url)
    if not r.ok:
        raise Exception(f"Request failed with status code {r.status_code}")

    content = r.content.decode()
    results = []
    for line in content.splitlines():
        if line.startswith("UP"):
            l = line.split("\t")
            if kingdom.lower() in l[3]:
                upcode = l[0]
                taxid = l[1]
                results.append( (kingdom, upcode, taxid) )
    return results


def retrieve(kingdom, upcode, taxid):
    url = f"{URL}/{kingdom}/{upcode}/{upcode}_{taxid}.fasta.gz"
    r = requests.get(url)
    if r.ok:
        destination = UP_DIR / kingdom / f"{upcode}.fasta.gz"
        with open(destination, "wb") as gzfile:
            gzfile.write(r.content)
    else:
        print(f"request failed: {url}")


def main(kingdom):
    inputs = get_list(kingdom)
    for params in inputs:
        retrieve(*params)


if __name__ == '__main__':
    main("Archaea")
