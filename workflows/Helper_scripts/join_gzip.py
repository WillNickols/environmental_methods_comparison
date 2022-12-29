import multiprocessing as mp
import glob
import os
import shutil
from pathlib import Path
import argparse
import subprocess

# python join_gzip.py --input input_folder --threads 20

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="path file")
parser.add_argument("--threads", help="threads", type=int)
args = parser.parse_args()

paths = Path("/" + (args.input).strip("/") + "/").rglob('*.fastq.gz')

files = [path.as_posix() for path in paths]

file_list = list(set([file.split("_R")[0] for file in files]))

def join_gzip(file):
    subprocess.run(['cat ' + file + '_R1.fastq.gz ' + file + '_R2.fastq.gz > ' + file + '.fastq.gz && rm ' + file + '_R1.fastq.gz && rm ' + file + '_R2.fastq.gz '], shell=True)

with mp.Pool(processes = args.threads) as p:
    p.map(join_gzip, file_list)
