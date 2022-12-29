#!/usr/bin/env python

from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import math
import shutil

workflow = Workflow(version="0.1", description="kraken workflow")
workflow.add_argument("database-folder", desc="Folder with the phylophlan database")
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
database_folder = args.database_folder

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
		if step == "kraken":
			time = 30 * n_gigabytes
		if step == "get_reads":
			time = 5 * n_gigabytes
	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "kraken":
			time = 15 * n_gigabytes
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
		if step == "kraken" or step == "get_reads":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
	if paired == "unpaired":
		if step == "kraken" or step == "get_reads":
			return [str(name + "." + input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "kraken":
		targets = output + name.split("/")[-1] + "_report_MPA_from_bracken.txt"
	if step == "get_reads":
		targets = output + name.split("/")[-1] + ".mapped_read_num.txt"
	return [str(targets)]

###########################
# run abundance profiling #
###########################

def kraken(name, paired):
	kraken_dir = scratch + "kraken/"
	if paired == "paired":
		command = '''{a} && {b} && {c}'''.format(
			a = "kraken2 --paired --db " + args.database_folder + " " + name + "_paired_1." + input_extension + " " + name + "_paired_2." + input_extension + " --threads " + str(cores) + " --output " + kraken_dir + name.split("/")[-1] + ".tsv --gzip-compressed --report " + output + name.split("/")[-1] + "_report_kraken.txt",
			b = "bash bracken -d " + args.database_folder + " -i " + output + name.split("/")[-1] + "_report_kraken.txt" + " -o " + output + name.split("/")[-1] + "_report_bracken.txt",
			c = "python kreport2mpa.py -r " + output + name.split("/")[-1] + "_report_kraken_bracken_species.txt -o " + output + name.split("/")[-1] + "_report_MPA_from_bracken.txt --display-header",
			)
	if paired == "unpaired":
		command = '''{a} && {b} && {c}'''.format(
			a = "kraken2 --db " + args.database_folder + " " + name + "." + input_extension + " --threads " + str(cores) + " --output " + kraken_dir + name.split("/")[-1] + ".tsv --gzip-compressed --report " + output + name.split("/")[-1] + "_report_kraken.txt",
			b = "bash bracken -d " + args.database_folder + " -i " + output + name.split("/")[-1] + "_report_kraken.txt" + " -o " + output + name.split("/")[-1] + "_report_bracken.txt",
			c = "python kreport2mpa.py -r " + output + name.split("/")[-1] + "_report_kraken_bracken_species.txt -o " + output + name.split("/")[-1] + "_report_MPA_from_bracken.txt --display-header",
			)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="kraken", paired=paired)[0]):
		workflow.add_task_gridable(actions=kraken(name, paired),
			depends=list_depends(name=name, step="kraken", paired=paired),
			targets=list_targets(name=name, step="kraken", paired=paired),
			time=calculate_time(name=name, step="kraken", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

combine_dep = [output + name.split("/")[-1] + "_report_MPA_from_bracken.txt" for name in names]
combine_targets = [output+"kraken_merged_tmp.tsv"]

workflow.add_task("python combine_mpa.py -i " + output + "*_report_MPA_from_bracken.txt -o [targets[0]]",
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
wrangling_dep.append(output+"kraken_merged_tmp.tsv")
wrangling_targets = [output+"kraken_merged.tsv"]

workflow.add_task("Rscript kraken_tax_wrangling.R -i " + output + " -o [targets[0]]",
	depends = wrangling_dep,
	targets = wrangling_targets
	)

################
# run workflow #
################

workflow.go()

#
