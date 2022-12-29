#!/usr/bin/env python

from anadama2 import Workflow
from pathlib import Path
import math
import os
import glob

workflow = Workflow(version="0.4", description="metaspades_test")

args = workflow.parse_args()

out = "/" + args.output.strip("/") + "/"
scratch = "/" + args.grid_scratch.strip("/") + "/"
partition = args.grid_partition

##############
# run GTDBTK #
##############

command = "python spades.py -1 /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/reSRR9320187_paired_1.fastq.gz -2 /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/reSRR9320187_paired_2.fastq.gz -s /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/reSRR9320187_unmatched_1.fastq.gz -s /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/reSRR9320187_unmatched_2.fastq.gz --meta -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/metaspades_test/ -t 12"

workflow.add_task_gridable(actions=command,
	time=60*24,
	mem=184000,
	cores=12,
	partition=partition
	)

workflow.go()
