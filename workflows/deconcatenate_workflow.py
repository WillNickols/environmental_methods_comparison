#!/usr/bin/env python

from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import math
import shutil
from pathlib import Path

workflow = Workflow(version="0.1", description="deconcatenate kneaddata files")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=1)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=3000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=1000)
args = workflow.parse_args()
input_extension = args.input_extension

# output
output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)

# scratch directory
scratch = "/" + args.grid_scratch.strip("/") + "/"

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time

# list the input fastq files
in_dir = args.input

paths = Path("/" + in_dir.strip("/") + "/").rglob('*.' + input_extension)

files = []
for path in paths:
	files.append(path.as_posix())

names = set(file.split('.' + input_extension)[0] for file in files)

#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step):
	n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024))
	if step == "deconcatenate":
		time = 20 * n_gigabytes
	return min([max_time, int(time)])

#################################
# function to list dependencies #
#################################

def list_depends(name, step):
    if step == "deconcatenate":
        return [str(name + "." + input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step):
	if step == "deconcatenate":
		return [str(output + name.split("/")[-1] + "_paired_1.fastq.gz"), str(output + name.split("/")[-1] + "_paired_2.fastq.gz"), str(output + name.split("/")[-1] + "_unmatched_1.fastq.gz"), str(output + name.split("/")[-1] + "_unmatched_2.fastq.gz")]

#############################
# function to run kneaddata #
#############################

def deconcatenate(name):
    if input_extension == "fastq.gz":
        unzipped_name = scratch + name.split("/")[-1] + ".fastq"
        command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {k}'''.format(
			a = "gunzip -c " + name + "." + input_extension + " > " + unzipped_name,
			b = "python deconcatenate.py " + unzipped_name + " " + scratch + name.split("/")[-1],
			c = "rm " + unzipped_name,
			d = "gzip " + scratch + name.split("/")[-1] + "_unmatched_1.fastq",
			e = "gzip " + scratch + name.split("/")[-1] + "_unmatched_2.fastq",
			f = "cat " + scratch + name.split("/")[-1] + "_paired_1.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > " + scratch + name.split("/")[-1] + "_1_sorted.fastq.gz",
			g = "cat " + scratch + name.split("/")[-1] + "_paired_2.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > " + scratch + name.split("/")[-1] + "_2_sorted.fastq.gz",
			h = "rm " + scratch + name.split("/")[-1] + "_paired_1.fastq",
			i = "rm " + scratch + name.split("/")[-1] + "_paired_2.fastq",
			j = "mv " + scratch + name.split("/")[-1] + "_1_sorted.fastq.gz " + scratch + name.split("/")[-1] + "_paired_1.fastq.gz",
			k = "mv " + scratch + name.split("/")[-1] + "_2_sorted.fastq.gz " + scratch + name.split("/")[-1] + "_paired_2.fastq.gz"
			)
    else:
        command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i}'''.format(
			a = "python deconcatenate.py " + name + "." + input_extension + " " + scratch + name.split("/")[-1],
			b = "gzip " + scratch + name.split("/")[-1] + "_unmatched_1.fastq",
			c = "gzip " + scratch + name.split("/")[-1] + "_unmatched_2.fastq",
			d = "cat " + scratch + name.split("/")[-1] + "_paired_1.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > " + scratch + name.split("/")[-1] + "_1_sorted.fastq.gz",
			e = "cat " + scratch + name.split("/")[-1] + "_paired_2.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip > " + scratch + name.split("/")[-1] + "_2_sorted.fastq.gz",
			f = "rm " + scratch + name.split("/")[-1] + "_paired_1.fastq",
			g = "rm " + scratch + name.split("/")[-1] + "_paired_2.fastq",
			h = "mv " + scratch + name.split("/")[-1] + "_1_sorted.fastq.gz " + scratch + name.split("/")[-1] + "_paired_1.fastq.gz",
			i = "mv " + scratch + name.split("/")[-1] + "_2_sorted.fastq.gz " + scratch + name.split("/")[-1] + "_paired_2.fastq.gz"
			)
    return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="deconcatenate")[0]):
		workflow.add_task_gridable(actions=deconcatenate(name),
			depends=list_depends(name=name, step="deconcatenate"),
			targets=list_targets(name=name, step="deconcatenate"),
			time=calculate_time(name=name, step="deconcatenate"),
			mem=memory,
			cores=cores,
			partition=partition
			)

####################
# run the workflow #
####################

workflow.go()

#
