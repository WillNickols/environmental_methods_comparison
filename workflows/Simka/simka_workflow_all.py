#!/usr/bin/env python

from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import math
import shutil

# Inputs must be gzipped
workflow = Workflow(version="0.1", description="simka workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=20000)
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("simka-path")
args = workflow.parse_args()

if args.simka_path.split("/")[-1] != "simka":
	raise ValueError("Simka path invalid")

# output
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

output_gunzip = output + "gunzip/"
if not os.path.isdir(output_gunzip):
	os.makedirs(output_gunzip)

# scratch directory
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"
scratch_in = scratch + "in/"
simka_scratch = scratch + "simka_tmp/"

if not os.path.isdir(scratch):
	os.makedirs(scratch)

if not os.path.isdir(scratch_in):
	os.makedirs(scratch_in)

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time

# list the input fastq files
with open(args.input) as f:
	inputs = f.readlines()

inputs_dict = {line.split(", ")[0]: line.split(", ")[1].strip("\n") for line in inputs}

i = 0
depends_list = []
all_names = []
all_paired = []
with open(output + "in.txt", "w") as f:
	for key, val in inputs_dict.items():
		if val == "paired" or val == "no_unmatched":
			paths = glob.glob(os.path.abspath(key.rstrip("/")) + "/*.fastq.gz")
			files = []
			for path in paths:
				files.append(path)
			names = set(file.split('_paired')[0].split('_unmatched')[0] for file in files)
			names = sorted(names)
			all_names.extend(names)
			all_paired.extend([val] * len(names))

			if val == "paired":
				for name in names:
					f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + "_paired_1.fastq; " + scratch_in + name.split("/")[-1] + "_paired_2.fastq; " + scratch_in + name.split("/")[-1] + "_unmatched_1.fastq; " + scratch_in + name.split("/")[-1] + "_unmatched_2.fastq\n")
					i += 1
			else:
				for name in names:
					f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + "_paired_1.fastq; " + scratch_in + name.split("/")[-1] + "_paired_2.fastq\n")
					i += 1

			tmp_list = [[str(scratch_in + name.split("/")[-1] + "_paired_1.fastq"), str(scratch_in + name.split("/")[-1] + "_paired_2.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_1.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_2.fastq")] for name in names]
			depends_list.extend([item for sublist in tmp_list for item in sublist])
				
		if val == "unpaired":
			paths = glob.glob(os.path.abspath(key.rstrip("/")) + "/*.fastq.gz")
			files = []
			for path in paths:
				files.append(path)
			names = set(file.split(".fastq.gz")[0] for file in files)
			names = sorted(names)

			all_names.extend(names)
			all_paired.extend([val] * len(names))

			for name in names:
				f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + ".fastq\n")
				i += 1

			depends_list.extend([str(scratch_in + name.split("/")[-1] + ".fastq") for name in names])

#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	time = 0
	if step == "simka":
		time = 15 + 15 * len(all_names)
	else:
		if paired == "paired" or paired == "no_unmatched":
			n_gigabytes = math.ceil(os.path.getsize(name + "_paired_1.fastq.gz") / (1024 * 1024 * 1024.0))
			if step == "gunzip":
				time = 20 * n_gigabytes
		if paired == "unpaired":
			n_gigabytes = math.ceil(os.path.getsize(name + ".fastq.gz") / (1024 * 1024 * 1024.0))
			if step == "gunzip":
				time = 10 * n_gigabytes
	if time > max_time:
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired" or paired == "no_unmatched":
		if step == "gunzip":
			return [str(name + "_paired_1.fastq.gz"), str(name + "_paired_2.fastq.gz"), str(name + "_unmatched_1.fastq.gz"), str(name + "_unmatched_2.fastq.gz")]
	if paired == "unpaired":
		if step == "gunzip":
			return [str(name + ".fastq.gz")]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "simka":
			return [output + "mat_abundance_braycurtis.csv.gz", output + "mat_abundance_jaccard.csv.gz"]
	if step == "gunzip":
		return [output_gunzip + name.split("/")[-1] + "_gunzip.done"]

##################
# copy and unzip #
##################

def gunzip(name, paired):
	if paired == "paired" or paired == "no_unmatched":
		command = '''{a} && {b} && {c}'''.format(
			a = "cp " + str(name + "_paired_1.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_paired_2.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_unmatched_1.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_unmatched_2.fastq.gz") + " " + scratch_in,
			b = "gunzip " + scratch_in + name.split("/")[-1] + "_paired_1.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_paired_2.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_unmatched_1.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_unmatched_2.fastq.gz",
			c = "touch " + output_gunzip + name.split("/")[-1] + "_gunzip.done"
			)
	if paired == "unpaired":
		command = '''{a} && {b} && {c}'''.format(
			a = "cp " + str(name + ".fastq.gz") + " " + scratch_in,
			b = "gunzip " + scratch_in + name.split("/")[-1] + ".fastq.gz",
			c = "touch " + output_gunzip + name.split("/")[-1] + "_gunzip.done"
			)
	return str(command)

for name, paired in zip(all_names, all_paired):
	if not os.path.isfile(list_targets(name=name, step="gunzip", paired=paired)[0]):
		workflow.add_task_gridable(actions=gunzip(name, paired),
			depends=list_depends(name=name, step="gunzip", paired=paired),
			targets=list_targets(name=name, step="gunzip", paired=paired),
			time=calculate_time(name=name, step="gunzip", paired=paired),
			mem=10000,
			cores=1,
			partition=partition
			)

#############
# run simka #
#############

if not os.path.isfile(list_targets(name="", step="simka", paired="")[0]):
	command = args.simka_path + " -in " + output + "in.txt -out " + scratch + " -out-tmp " + simka_scratch + " -max-memory " + str(memory) + " -nb-cores " + str(cores)
	workflow.add_task_gridable(actions=command,
		depends=[output_gunzip + name.split("/")[-1] + "_gunzip.done" for name in all_names],
		targets=list_targets(name="", step="simka", paired=""),
		time=calculate_time(name="", step="simka", paired=""),
		mem=memory,
		cores=cores,
		partition=partition
		)

##########################
# delete temporary files #
##########################

workflow.add_task(actions="rm -r " + scratch + " && gunzip " + output + "mat_abundance_braycurtis.csv.gz && gunzip " + output + "mat_abundance_jaccard.csv.gz",
	depends=list_targets(name="", step="simka", paired=""),
	targets=[output + "mat_abundance_braycurtis.csv", output + "mat_abundance_jaccard.csv"]
	)

################
# run workflow #
################

workflow.go()
