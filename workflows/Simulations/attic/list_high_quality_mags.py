# To run: $ python list_high_quality_mags.py --qa checkm_qa_and_n50.tsv --bins-path bins/ --output out.tsv

import multiprocessing as mp
import pandas as pd
import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--qa", help="checkm_qa_and_n50.tsv")
parser.add_argument("--bins-path", help="bins/")
parser.add_argument("--output", help="out.tsv")
args = parser.parse_args()

df = pd.read_csv(args.qa, sep="\t")

bin_list = df.loc[df["keep"] == "keep", "bin_id"].tolist()
bin_list = [args.bins_path + bin.split(".bin.")[0] + "/bins/" + bin + ".fa" for bin in bin_list]

pd.DataFrame(bin_list).to_csv(args.output, mode='a', sep="\t", header=False, index=False)
