#!/usr/bin/env python

import pickle
import bz2
from anadama2 import Workflow
import os
import csv

workflow = Workflow(version="0.1", description="Get all MPA2 taxonomic groups")
args = workflow.parse_args()

# output
output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)

db2 = pickle.load(bz2.BZ2File(args.input, 'r'))
print("Loaded")

taxa = set(db2["taxonomy"].keys())

with open(output+'MPA2Taxonomy.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerows([item] for item in list(taxa))

print("Done")
