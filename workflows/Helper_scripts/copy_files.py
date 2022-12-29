# To run: $ python copy_files.py --paths $PATHS_DOC --threads 4 --output $OUTPUT_FOLDER

import multiprocessing as mp
import pandas as pd
import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--paths", help="path file")
parser.add_argument("--threads", help="threads", type=int)
parser.add_argument("--output", help="output")
args = parser.parse_args()

output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
    os.makedirs(output)

df = pd.read_csv(args.paths, header = None)

in_files = (df.iloc[:, 0]).tolist()
out_files = [output + s.split("/")[-1] for s in in_files]

parallel_list = [(in_file, out_file) for in_file, out_file in zip(in_files, out_files)]

def copy_file(in_file, out_file):
    subprocess.run(['cp ' + in_file + ' ' + out_file], shell=True)

with mp.Pool(processes = args.threads) as p:
    p.starmap(copy_file, parallel_list)
