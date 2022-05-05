from pathlib import Path

ROOT = Path(__file__).parent.parent


# sequence databases
DATA_DIR = ROOT / 'data'
UP_DIR = DATA_DIR / 'uniprot'
NCBI_DIR = DATA_DIR / 'ncbi'

# search queries and results
SEARCH_DIR = ROOT / "searches"
QUERY_DIR = SEARCH_DIR / 'queries'
HALO_DIR = SEARCH_DIR / "halocyanin"
CYTO_DIR = SEARCH_DIR / "cytochrome"

# target fasta files
TARGETS = ROOT / "targets"

# average proteome data
PROTEOMEDATA = ROOT / "analysis" / "proteomedata"

# prosite results
PROSITES = ROOT / "prosites"
