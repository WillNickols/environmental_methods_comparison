import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import os
import glob
import multiprocessing as mp
import argparse
from pathlib import Path

# Inputs needed:
# new_sgbs/sgbs/sgbs/ path
# simulations/genomes_dir/ path
# profile.tsv which mediates all BIN labels

parser = argparse.ArgumentParser()
parser.add_argument("--sgb-path", help="sgbs/sgbs/ path")
parser.add_argument("--genomes-dir", help="genomes_dir/ path")
parser.add_argument("--profile", help="profile.tsv")
parser.add_argument("--out-file", help="out.tsv")
parser.add_argument("--threads", help="threads to use", type=int, default = 1)
args = parser.parse_args()

sgb_path = args.sgb_path
genomes_dir = args.genomes_dir
out_file = args.out_file
threads = args.threads

paths = Path(os.path.abspath(sgb_path.rstrip("/")) + "/").rglob("*.fa")

files = []
for path in paths:
    files.append(str(path))

paths = Path(os.path.abspath(genomes_dir.rstrip("/")) + "/known/genomes/").rglob("*.fa")

for path in paths:
    files.append(str(path))

paths = Path(os.path.abspath(genomes_dir.rstrip("/")) + "/new/genomes/").rglob("*.fa")

for path in paths:
    files.append(str(path))

def calculate_gc(fasta):
    total_length = 0
    for record in SeqIO.parse(fasta,'fasta'):
        total_length += len(record)

    return (fasta.split("/")[-1].split(".fa")[0], total_length)

with mp.Pool(processes = threads) as p:
    gc_out = p.map(calculate_gc, files)

os.makedirs(out_file.rsplit("/", 1)[0], exist_ok=True)

length_df = pd.DataFrame(gc_out, columns=['name', 'length'])

profile = pd.read_csv(args.profile, sep="\t")
profile['name'] = profile['filepath'].str.split("/").str[-1].str.split(".fa").str[0]

length_df = length_df.merge(profile, on='name', how='left')
length_df = length_df.drop(columns=['filepath', 'NCBI_ID', 'novelty_category'])

length_df.to_csv(out_file, sep='\t', index=False)