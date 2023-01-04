#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil

workflow = Workflow(version="0.1", description="mOTUs3 workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=45000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument(name="paired", default="paired")
args = workflow.parse_args()
input_extension = args.input_extension

# output
output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

# scratch directory
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time
paired = args.paired

# list the input fastq files
in_dir = args.input

if paired == "paired":
	paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/" + '*.' + input_extension)
	files = []
	for path in paths:
		files.append(path)
	names = set(file.split('_paired')[0].split('_unmatched')[0] for file in files)

if paired == "unpaired":
	paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/" + '*.' + input_extension)
	files = []
	for path in paths:
		files.append(path)
	names = set(file.split("." + input_extension)[0] for file in files)


#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	time = 0
	if paired == "paired":
		n_gigabytes = math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "mOTUs":
			time = 20 * n_gigabytes
	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "mOTUs":
			time = 12 * n_gigabytes
	if input_extension == "fastq.gz":
		time = time * 2
	if time > max_time:
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired":
		if step == "mOTUs":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
	if paired == "unpaired":
		if step == "mOTUs":
			return [str(name + "." + input_extension)]
	return [str(depends)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "mOTUs":
		target = output + "main/" + name.split("/")[-1] + ".txt"
	return [str(target)]


###########################
# function to call mOTUs3 #
###########################

def mOTUs(name, paired):
	if paired == "paired":
		command = '''{a}'''.format(
			a = "motus profile -f " + name + "_paired_1." + input_extension + " -r " + name + "_paired_2." + input_extension + " -t " + str(cores) + " -A -o [targets[0]] -n " + name.split("/")[-1]
			)
	if paired == "unpaired":
		command = '''{a}'''.format(
			a = "motus profile -s " + name + "." + input_extension + " -t " + str(cores) + " -A -o [targets[0]] -n " + name.split("/")[-1]
			)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="mOTUs", paired=paired)[0]):
		workflow.add_task_gridable(actions=mOTUs(name, paired),
			depends=list_depends(name=name, step="mOTUs", paired=paired),
			targets=list_targets(name=name, step="mOTUs", paired=paired),
			time=calculate_time(name=name, step="mOTUs", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

#############################
# function to run MetaBAT 2 #
#############################

def merge():
	if not os.path.isdir(output + "merged/"):
		os.makedirs(output + "merged/")

	command = '''{a}'''.format(
		a = "motus merge -d " + output + "main/ -o " + output + "merged/all_sample_profiles.txt"
		)

	return str(command)

depends = [list_targets(name=name, step="mOTUs", paired=paired)[0] for name in names]

if not os.path.isfile(output + "merged/all_sample_profiles.txt"):
	workflow.add_task_gridable(actions=merge(),
		depends=depends,
		targets=output + "merged/all_sample_profiles.txt",
		time=5,
		mem=memory,
		cores=1,
		partition=partition
		)

wrangling_dep = [output + "merged/all_sample_profiles.txt"]
wrangling_targets = [output + "mOTUs3_merged.tsv"]

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

workflow.add_task("Rscript " + this_folder + "mOTUs3_tax_wrangling.R -i [depends[0]] -o [targets[0]]",
	depends = wrangling_dep,
	targets = wrangling_targets
	)

####################
# run the workflow #
####################

workflow.go()

#
