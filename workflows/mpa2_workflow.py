#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil

workflow = Workflow(version="0.1", description="MPA 2 workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=45000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument(name="paired", default="paired")
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
paired = args.paired

# list the input fastq files
in_dir = args.input

if paired == "paired":
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)

	files = []
	for path in paths:
		files.append(path)

	names = set(file.split('_paired')[0].split('_unmatched')[0] for file in files)

if paired == "unpaired":
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)

	files = []
	for path in paths:
		files.append(path)

	names = set(file.split("." + input_extension)[0] for file in files)


#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	time = 1	

	if paired == "paired":
		n_gigabytes = math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "metaphlan2":
			time = 40 * n_gigabytes

	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "metaphlan2":
			time = 20 * n_gigabytes

	if time > max_time:
		time = max_time

	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired":
		if step == "metaphlan2":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension), str(name + "_unmatched_1." + input_extension), str(name + "_unmatched_2." + input_extension)]
	if paired == "unpaired":
		if step == "metaphlan2":
			return [str(name+"."+input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if paired == "paired":
		if step == "metaphlan2":
			return [str(output + name.split("/")[-1] + "_taxonomic_profile.tsv")]
	if paired == "unpaired":
		if step == "metaphlan2":
			return [str(output + name.split("/")[-1] + "_taxonomic_profile.tsv")]

##############################
# function to run metaphlan2 #
##############################

def metaphlan2(name, paired):
	if paired == "paired":
		if input_extension == "fastq.gz":
			command = '''{a}'''.format(
				a = "metaphlan2.py <(zcat " + name + "_paired_1." + input_extension + " " + name + "_paired_2." + input_extension + " " + name + "_unmatched_1." + input_extension + " " + name + "_unmatched_2." + input_extension + ") --input_type fastq --output_file " + scratch + name.split("/")[-1] + "_taxonomic_profile.tsv --nproc " + str(cores) + " --no_map --tmp_dir " + scratch
				)
	if paired == "unpaired":
		if input_extension == "fastq.gz":
			command = '''{a}'''.format(
				a = "metaphlan2.py <(zcat " + name + "." + input_extension + ") --input_type fastq --output_file " + scratch + name.split("/")[-1] + "_taxonomic_profile.tsv --nproc " + str(cores) + " --no_map --tmp_dir " + scratch
				)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="metaphlan2", paired=paired)[0]):
		workflow.add_task_gridable(actions=metaphlan2(name, paired=paired),
			depends=list_depends(name=name, step="metaphlan2", paired=paired),
			targets=list_targets(name=name, step="metaphlan2", paired=paired),
			time=calculate_time(name=name, step="metaphlan2", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

merge_depends = []
for name in names:
	merge_depends.append(list_targets(name=name, step="metaphlan2", paired=paired)[0])

final_output = [output + "/merged/metaphlan2_taxonomic_profiles.tsv"]

command = '''{a}'''.format(
	a = "merge_metaphlan_tables.py " + output + "*.tsv > [targets[0]]"
	)

workflow.add_task_gridable(actions=command,
	depends=merge_depends,
	targets=final_output,
	time=5,
	mem=memory,
	cores=1,
	partition=partition
	)

####################
# run the workflow #
####################

workflow.go()

#
