#!/usr/bin/env python

from anadama2 import Workflow
from pathlib import Path
import os
import math
import shutil

workflow = Workflow(version="0.1", description="Kneaddata workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=45000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument(name="paired", default="paired")
workflow.add_argument(name="reference-db", default="/n/huttenhower_lab/data/kneaddata_databases/")
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
	paths = Path(os.path.abspath(in_dir.rstrip("/")) + "/").rglob('*.' + input_extension)

	files = []
	for path in paths:
		files.append(path.as_posix())

	names = set(file.rsplit('_1', 1)[0] for file in files if "_1." + input_extension in file)

if paired == "unpaired":
	paths = Path(os.path.abspath(in_dir.rstrip("/")) + "/").rglob('*.' + input_extension)

	files = []
	for path in paths:
		files.append(path.as_posix())

	names = set(file.split("." + input_extension)[0] for file in files)


#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	print(name)
	if paired == "paired":
		n_gigabytes = math.ceil(os.path.getsize(name + "_1." + input_extension) / (1024 * 1024 * 1024))
		if step == "kneaddata":
			time = 80 * n_gigabytes
		return int(time)

	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024))
		if step == "kneaddata":
			time = 30 * n_gigabytes
		return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired":
		if step == "kneaddata":
			return [str(name+"_1."+input_extension), str(name+"_2."+input_extension)]
	if paired == "unpaired":
		if step == "kneaddata":
			return [str(name+"."+input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if paired == "paired":
		if step == "kneaddata":
			return [str(output + name.split("/")[-1] + "_paired_1.fastq.gz"), str(output + name.split("/")[-1] + "_paired_2.fastq.gz"), str(output + name.split("/")[-1] + "_unmatched_1.fastq.gz"), str(output + name.split("/")[-1] + "_unmatched_2.fastq.gz")]
	if paired == "unpaired":
		if step == "kneaddata":
			targets0 = output + name.split("/")[-1] + ".fastq.gz"
			return [str(targets0)]

#############################
# function to run kneaddata #
#############################

def kneaddata(name, paired):
	if paired == "paired":
		command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i}'''.format(
			a = "kneaddata --input " + name + "_1." + input_extension + " --output " + scratch + " --threads " + str(cores) + " --output-prefix " + name.split("/")[-1] + " --input " + name + "_2." + input_extension + " --reference-db " + args.reference_db + " --serial --run-trf  --remove-intermediate-output",
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
	if paired == "unpaired":
		command = '''{a} && {b}'''.format(
			a = "kneaddata --input " + name + "." + input_extension + " --output " + scratch + " --threads " + str(cores) + " --output-prefix " + name.split("/")[-1] + " --reference-db " + args.reference_db + " --serial --run-trf  --remove-intermediate-output",
			b = "gzip " + scratch + name.split("/")[-1] + ".fastq"
			)
		return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="kneaddata", paired=paired)[0]):
		workflow.add_task_gridable(actions=kneaddata(name, paired=paired),
			depends=list_depends(name=name, step="kneaddata", paired=paired),
			targets=list_targets(name=name, step="kneaddata", paired=paired),
			time=calculate_time(name=name, step="kneaddata", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

####################
# run the workflow #
####################

workflow.go()

#
