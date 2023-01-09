#!/usr/bin/env python

from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import math
import shutil

try:
	centrifuge_path = os.environ["CENTRIFUGE_PATH"].rstrip("/") + "/"
except:
	raise ValueError("Centrifuge path is not set (CENTRIFUGE_PATH)")

workflow = Workflow(version="0.1", description="centrifuge workflow")
workflow.add_argument("x", desc="-x parameter for centrifuge", default = centrifuge_path + "database/abv")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=8)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=40000)
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

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

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
		if step == "centrifuge":
			time = 10 + 16 * n_gigabytes
		if step == "get_reads":
			time = 5 * n_gigabytes
	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "centrifuge":
			time = 10 + 20 * n_gigabytes
		if step == "get_reads":
			time = 3 * n_gigabytes
	if time > max_time:
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired":
		if step == "centrifuge" or step == "get_reads":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
	if paired == "unpaired":
		if step == "centrifuge" or step == "get_reads":
			return [str(name + "." + input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "centrifuge":
		targets = output + name.split("/")[-1] + "_report_MPA.txt"
	if step == "get_reads":
		targets = output + name.split("/")[-1] + ".mapped_read_num.txt"
	return [str(targets)]

###########################
# run abundance profiling #
###########################

def centrifuge(name, paired):
	if paired == "paired":
		command = '''{a} && {b} && {c} && {d}'''.format(
			a = centrifuge_path + "centrifuge -x " + str(args.x) + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -U " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " --report-file " + output + name.split("/")[-1] + "_centrifuge_report.tsv " + " -S " + output + name.split("/")[-1] + "_individual_classification.txt" + " -p " + str(cores),
			b = centrifuge_path + "centrifuge-kreport" + " -x " + str(args.x) + " " + output + name.split("/")[-1] + "_individual_classification.txt > " + output + name.split("/")[-1] + "_report_centrifuge_species.txt",
			c = "python " + centrifuge_path + "kreport2mpa.py -r " + output + name.split("/")[-1] + "_report_centrifuge_species.txt -o " + output + name.split("/")[-1] + "_report_MPA.txt --display-header",
			d = "rm " + output + name.split("/")[-1] + "_centrifuge_report.tsv"
			)
	if paired == "unpaired":
		command = '''{a} && {b} && {c} && {d}'''.format(
			a = centrifuge_path + "centrifuge -x " + str(args.x) + " -U " + name + "." + input_extension + " --report-file " + output + name.split("/")[-1] + "_centrifuge_report.tsv " + " -S " + output + name.split("/")[-1] + "_individual_classification.txt" + " -p " + str(cores),
			b = centrifuge_path + "centrifuge-kreport" + " -x " + str(args.x) + " " + output + name.split("/")[-1] + "_individual_classification.txt > " + output + name.split("/")[-1] + "_report_centrifuge_species.txt",
			c = "python " + centrifuge_path + "kreport2mpa.py -r " + output + name.split("/")[-1] + "_report_centrifuge_species.txt -o " + output + name.split("/")[-1] + "_report_MPA.txt --display-header",
			d = "rm " + output + name.split("/")[-1] + "_centrifuge_report.tsv"
			)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="centrifuge", paired=paired)[0]):
		workflow.add_task_gridable(actions=centrifuge(name, paired),
			depends=list_depends(name=name, step="centrifuge", paired=paired),
			targets=list_targets(name=name, step="centrifuge", paired=paired),
			time=calculate_time(name=name, step="centrifuge", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

combine_dep = [output + name.split("/")[-1] + "_report_MPA.txt" for name in names]
combine_targets = [output+"centrifuge_merged_tmp.tsv"]

workflow.add_task("python " + centrifuge_path + "combine_mpa.py -i " + output + "*_report_MPA.txt -o [targets[0]]",
	depends = combine_dep,
	targets = combine_targets
	)

def get_reads(name, paired):
	if paired == "paired":
		if input_extension == "fastq.gz":
			command = '''{a} && {b} && {c} && {d}'''.format(
				a = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				b = "unpaired1=$(echo $(zcat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				c = "unpaired2=$(echo $(zcat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				d = "echo $((paired1+unpaired1+unpaired2)) > [targets[0]]"
				)

	if paired == "unpaired":
		if input_extension == "fastq.gz":
			command = '''{a}'''.format(
				a = "echo $(zcat " + name + "." + input_extension + "|wc -l)/4|bc > [targets[0]]"
				)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="get_reads", paired=paired)[0]):
		workflow.add_task_gridable(actions=get_reads(name, paired),
			depends=list_depends(name=name, step="get_reads", paired=paired),
			targets=list_targets(name=name, step="get_reads", paired=paired),
			time=calculate_time(name=name, step="get_reads", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

wrangling_dep = [list_targets(name=name, step="get_reads", paired=paired)[0] for name in names]
wrangling_dep.append(output+"centrifuge_merged_tmp.tsv")
wrangling_targets = [output+"centrifuge_merged.tsv"]

workflow.add_task("Rscript " + this_folder + "centrifuge_tax_wrangling.R -i " + output + " -o [targets[0]]",
	depends = wrangling_dep,
	targets = wrangling_targets
	)

################
# run workflow #
################

workflow.go()

#
