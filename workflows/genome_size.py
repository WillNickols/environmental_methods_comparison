#!/usr/bin/env python

# Use: python genome_size.py sgb_directory output_path

from pathlib import Path
import pandas as pd
import sys as sys

sgb_directory = "/" + str(sys.argv[1]).strip("/") + "/"
output_path = "/" + str(sys.argv[2]).strip("/") + "/"

paths = Path(sgb_directory).rglob("*.fa")

genomes = []

for path in paths:
	genomes.append(path.as_posix())

sizes = []

for genome in genomes:
	fasta = open(genome, "r")
	n = 0
	with open(genome) as fasta:
		for sequence in fasta:
			if ">" not in sequence:
				n += len(sequence)
	sizes.append(pd.DataFrame(data = {"sgb": [genome], "bp": [n]}))

df = pd.concat(sizes)

df.to_csv(output_path + "sgb_sizes.tsv", sep="\t", index=False)
