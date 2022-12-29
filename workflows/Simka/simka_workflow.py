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
workflow.add_argument(name="paired", default="paired")
workflow.add_argument("no-unmatched")
workflow.add_argument("simka-path")
args = workflow.parse_args()

if args.simka_path.split("/")[-1] != "simka":
	raise ValueError("Simka path invalid")

# output
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

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
paired = args.paired

# list the input fastq files
in_dir = args.input

if paired == "paired":
	paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/*.fastq.gz")
	files = []
	for path in paths:
		files.append(path)
	names = set(file.split('_paired')[0].split('_unmatched')[0] for file in files)

if paired == "unpaired":
	paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/*.fastq.gz")
	files = []
	for path in paths:
		files.append(path)
	names = set(file.split(".fastq.gz")[0] for file in files)

names = sorted(names)

#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	time = 0
	if step == "simka":
		if paired == "paired":
			time = 15 + 10 * len(names) * math.ceil(os.path.getsize(name + "_paired_1.fastq.gz") / (1024 * 1024 * 1024.0))
		else:
			time = 15 + 5 * len(names) * math.ceil(os.path.getsize(name + ".fastq.gz") / (1024 * 1024 * 1024.0))
	else:
		if paired == "paired":
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
	if paired == "paired":
		if step == "gunzip":
			return [str(name + "_paired_1.fastq.gz"), str(name + "_paired_2.fastq.gz"), str(name + "_unmatched_1.fastq.gz"), str(name + "_unmatched_2.fastq.gz")]
		if step == "simka":
			depends_list = [[str(scratch_in + name.split("/")[-1] + "_paired_1.fastq"), str(scratch_in + name.split("/")[-1] + "_paired_2.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_1.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_2.fastq")] for name in names]
			return [item for sublist in depends_list for item in sublist]
	if paired == "unpaired":
		if step == "gunzip":
			return [str(name + ".fastq.gz")]
		if step == "simka":
			return [str(scratch_in + name.split("/")[-1] + ".fastq") for name in names]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if paired == "paired":
		if step == "gunzip":
			return [str(scratch_in + name.split("/")[-1] + "_paired_1.fastq"), str(scratch_in + name.split("/")[-1] + "_paired_2.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_1.fastq"), str(scratch_in + name.split("/")[-1] + "_unmatched_2.fastq")]
		if step == "simka":
			return [output + "mat_abundance_braycurtis.csv.gz", output + "mat_abundance_jaccard.csv.gz"]
	if paired == "unpaired":
		if step == "gunzip":
			return [str(scratch_in + name.split("/")[-1] + ".fastq")]
		if step == "simka":
			return [output + "mat_abundance_braycurtis.csv.gz", output + "mat_abundance_jaccard.csv.gz"]

##################
# copy and unzip #
##################

def gunzip(name, paired):
	if paired == "paired":
		command = '''{a} && {b}'''.format(
			a = "cp " + str(name + "_paired_1.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_paired_2.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_unmatched_1.fastq.gz") + " " + scratch_in + " && cp " + str(name + "_unmatched_2.fastq.gz") + " " + scratch_in,
			b = "gunzip " + scratch_in + name.split("/")[-1] + "_paired_1.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_paired_2.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_unmatched_1.fastq.gz && gunzip " + scratch_in + name.split("/")[-1] + "_unmatched_2.fastq.gz"
			)
	if paired == "unpaired":
		command = '''{a} && {b}'''.format(
			a = "cp " + str(name + ".fastq.gz") + " " + scratch_in,
			b = "gunzip " + scratch_in + name.split("/")[-1] + ".fastq.gz"
			)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="gunzip", paired=paired)[0]):
		workflow.add_task(actions=gunzip(name, paired),
			depends=list_depends(name=name, step="gunzip", paired=paired),
			targets=list_targets(name=name, step="gunzip", paired=paired)
			)

#####################
# write in.txt file #
#####################

if not args.no_unmatched:
	with open(scratch + "in.txt", "w") as f:
		if paired == "paired":
			for i, name in enumerate(names):
				f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + "_paired_1.fastq; " + scratch_in + name.split("/")[-1] + "_paired_2.fastq; " + scratch_in + name.split("/")[-1] + "_unmatched_1.fastq; " + scratch_in + name.split("/")[-1] + "_unmatched_2.fastq\n")
		else:
			for i, name in enumerate(names):
				f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + ".fastq\n")
else:
	with open(scratch + "in.txt", "w") as f:
		if paired == "paired":
			for i, name in enumerate(names):
				f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + "_paired_1.fastq; " + scratch_in + name.split("/")[-1] + "_paired_2.fastq\n")
		else:
			for i, name in enumerate(names):
				f.write("ID" + str(i) + ": " + scratch_in + name.split("/")[-1] + ".fastq\n")

#############
# run simka #
#############

if not os.path.isfile(list_targets(name="", step="simka", paired=paired)[0]):
	command = args.simka_path + " -in " + scratch + "in.txt -out " + scratch + " -out-tmp " + simka_scratch + " -max-memory " + str(memory) + " -nb-cores " + str(cores)
	workflow.add_task_gridable(actions=command,
		depends=list_depends(name=name, step="simka", paired=paired),
		targets=list_targets(name=name, step="simka", paired=paired),
		time=calculate_time(name=name, step="simka", paired=paired),
		mem=memory,
		cores=cores,
		partition=partition
		)

##########################
# delete temporary files #
##########################

workflow.add_task(actions="rm -r " + scratch + " && gunzip " + output + "mat_abundance_braycurtis.csv.gz && gunzip " + output + "mat_abundance_jaccard.csv.gz",
	depends=list_targets(name=name, step="simka", paired=paired),
	targets=[output + "mat_abundance_braycurtis.csv", output + "mat_abundance_jaccard.csv"]
	)

################
# run workflow #
################

workflow.go()
