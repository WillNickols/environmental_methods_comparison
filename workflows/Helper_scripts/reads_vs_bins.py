#!/usr/bin/env python

from anadama2 import Workflow
from pathlib import Path
import os
import shutil
from collections import Counter
import pandas as pd

workflow = Workflow(version="0.1", description="Reads vs bins")
workflow.add_argument("checkm-dir")
args = workflow.parse_args()

checkm_directory = args.checkm_dir
output = args.output

paths = Path(os.path.abspath(checkm_directory.rstrip("/")) + "/").rglob('*.bin.*')

files = []
for path in paths:
	files.append(path.as_posix())

files = [file for file in files if not "lowDepth" in file and not "tooShort" in file]

names = [file.split("/")[-1].split(".bin.")[0] for file in files]
counts = Counter(names)
df = pd.DataFrame(list(counts.items()), columns=['ID', 'bins'])
df['reads'] = [None] * len(df.index)

paths = Path(os.path.abspath(checkm_directory.rstrip("/")) + "/").rglob('*.mapped_read_num.*')

files = []
for path in paths:
	files.append(path.as_posix())

df2 = pd.DataFrame(columns=['ID', 'reads'])

for file in files:
	name = (file.split(".")[0]).split("/")[-1]
	with open(file, 'r') as file_open:
		file_open.readline()
		tmp = file_open.readline()
		df2.loc[len(df2.index)] = [name, int(tmp)]
		df.loc[df.ID == name, 'reads'] = int(tmp)

df['bins'] = df['bins'].fillna(0)

# output
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)
df.to_csv(output+"reads_vs_bins.tsv", sep='\t', index=False)
df2.to_csv(output+"reads.tsv", sep='\t', index=False)
