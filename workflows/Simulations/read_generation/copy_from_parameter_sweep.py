import argparse
import multiprocessing as mp
import os
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--output", "-o", help="output")
parser.add_argument("--threads", help="threads",type=int)
args = parser.parse_args()

output = args.output
threads = args.threads

in_dir = output + "original/"
out_dir = output + "reorganized/"
if not os.path.isdir(out_dir + "reads/"):
    os.makedirs(out_dir + "reads/")

if not os.path.isdir(out_dir + "tax/"):
    os.makedirs(out_dir + "tax/")

paths = Path(os.path.abspath(in_dir.rstrip("/")) + "/").rglob("*anonymous_reads.fq.gz")

files = [str(path) for path in paths]

def move_fastq(fastq):
    fastq_pieces = fastq.split("/")
    new_file = out_dir + "reads/" + fastq_pieces[-7] + ".fastq.gz"
    os.rename(fastq, new_file)

with mp.Pool(processes = threads) as p:
    p.map(move_fastq, files)

paths = Path(os.path.abspath(in_dir.rstrip("/")) + "/").rglob("*profile_*.txt")

files = [str(path) for path in paths]

def move_tax(tax):
    tax_pieces = tax.split("/")
    new_file = out_dir + "tax/" + tax_pieces[-5] + ".txt"
    os.rename(tax, new_file)

with mp.Pool(processes = threads) as p:
    p.map(move_tax, files)
