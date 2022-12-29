import sys
from Bio import SeqIO
from Bio import Seq
import pandas as pd
import os
import glob
import multiprocessing as mp
import argparse
from pathlib import Path
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument("--in-dir", help="input directory")
parser.add_argument("--out-dir", help="output directory")
parser.add_argument("--rate", help="mutation frequency", type=float, default=0.001)
parser.add_argument("--threads", help="threads to use", type=int, default = 1)
args = parser.parse_args()

in_dir = args.in_dir
out_dir = args.out_dir
threads = args.threads
rate = args.rate

out_dir = os.path.abspath(out_dir.rstrip("/")) + "/"
if not os.path.isdir(out_dir):
	os.makedirs(out_dir)

paths = Path(os.path.abspath(in_dir.rstrip("/")) + "/").rglob("*.fa")
files = [str(path) for path in paths]

def mutate(fasta):
    sequences=[]
    for record in SeqIO.parse(fasta,'fasta'):
        indexes = np.random.randint(low = 0, high = len(record.seq), size = int(len(record.seq) * rate))
        sequence = list(record.seq)
        for i in indexes:
            sequence[i] = random.choice([x for x in "ACTG" if x != sequence[i].upper()])
        record.seq = Seq.Seq(''.join(sequence))
        sequences.append(record)
    SeqIO.write(sequences,out_dir + fasta.split("/")[-1],"fasta")

with mp.Pool(processes = threads) as p:
    gc_out = p.map(mutate, files)

with open(out_dir + "mutate.done", 'w') as fp:
    pass
