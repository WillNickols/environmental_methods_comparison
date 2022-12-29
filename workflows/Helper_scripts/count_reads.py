import subprocess
import multiprocessing as mp
import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--in-dir", help="input directory")
parser.add_argument("--threads", help="threads to train on", default = 4, type=int)
parser.add_argument("--output", help="output")
args = parser.parse_args()

fastqs = [str(f) for f in Path(args.in_dir).glob('*.fastq.gz')]

def count_reads(file):
    return (subprocess.run(['echo $(zcat ' + file + '|wc -l)/4|bc'], stdout=subprocess.PIPE, shell=True).stdout.decode("utf-8").split("\n")[0])

with mp.Pool(processes = args.threads) as p:
    read_nums = p.map(count_reads, fastqs)

pd.DataFrame.from_dict({"sample": fastqs, "reads": read_nums}).to_csv(args.output, sep='\t')
