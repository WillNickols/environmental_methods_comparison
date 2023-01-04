#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil

workflow = Workflow(version="0.1", description="metaxa2 workflow")
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
	time = 1

	if paired == "paired":
		n_gigabytes = math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "metaxa":
			time = 100 * n_gigabytes

	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "metaxa":
			time = 50 * n_gigabytes

	if step == "rewrite":
		time = 5

	if step == "collect":
		time = 5

	if input_extension == "fastq.gz":
		time = time * 5

	if time > max_time:
		time = max_time

	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if paired == "paired":
		if step == "metaxa":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension), str(name + "_unmatched_1." + input_extension), str(name + "_unmatched_2." + input_extension)]
	if paired == "unpaired":
		if step == "metaxa":
			return [str(name + "." + input_extension)]
	if step == "rewrite":
		depends = output + name.split("/")[-1] + ".taxonomy.txt"
	if step == "collect":
		depends = output + name.split("/")[-1] + ".level_1.txt"
	return [str(depends)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "metaxa":
		targets = output + name.split("/")[-1] + ".taxonomy.txt"
	if step == "rewrite":
		targets = output + name.split("/")[-1] + ".level_1.txt"
	if step == "collect":
		targets = output + "merged.tsv"
	return [str(targets)]

###########################
# function to call metaxa #
###########################

def metaxa(name, paired):
	if not os.path.isdir(output + name.split("/")[-1] + "/"):
		os.makedirs(output + name.split("/")[-1] + "/")
	if paired == "paired":
		if input_extension == "fastq.gz":
			command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			    a = "mkdir -p " + scratch + name.split("/")[-1] + "/",
				b = "zcat " + name + "_paired_1." + input_extension + " > " + scratch + name.split("/")[-1] + "_paired_1.fastq",
				c = "zcat " + name + "_paired_2." + input_extension + " > " + scratch + name.split("/")[-1] + "_paired_2.fastq",
				d = "metaxa2 -1 " + scratch + name.split("/")[-1] + "_paired_1.fastq -2 " + scratch + name.split("/")[-1] + "_paired_2.fastq -o " + scratch + name.split("/")[-1] + " --paired-fasta --cpu " + str(cores) + " --temp " + scratch[:-1],
				e = "rm " + scratch + name.split("/")[-1] + "_paired_1.fastq",
				f = "rm " + scratch + name.split("/")[-1] + "_paired_2.fastq",
				)
	if paired == "unpaired":
		if input_extension == "fastq.gz":
			command = '''{a} && {b} && {c}'''.format(
				a = "zcat " + name + "." + input_extension + " | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > " + scratch + name.split("/")[-1] + ".fastq",
			    b = "mkdir -p " + scratch + name.split("/")[-1] + "/",
				c = "metaxa2 -i " + scratch + name.split("/")[-1] + ".fastq -o " + scratch + name.split("/")[-1] + " --cpu " + str(cores)
				)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="metaxa", paired=paired)[0]):
		workflow.add_task_gridable(actions=metaxa(name=name, paired=paired),
			depends=list_depends(name=name, step="metaxa", paired=paired),
			targets=list_targets(name=name, step="metaxa", paired=paired),
			time=calculate_time(name=name, step="metaxa", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

################################
# function to rewrite taxonomy #
################################

def rewrite(name):
	command = '''{a}'''.format(
		a = "metaxa2_ttt -i [depends[0]] -o " + scratch + name.split("/")[-1]
		)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="rewrite", paired=paired)[0]):
		workflow.add_task_gridable(actions=rewrite(name),
			depends=list_depends(name=name, step="rewrite", paired=paired),
			targets=list_targets(name=name, step="rewrite", paired=paired),
			time=calculate_time(name=name, step="rewrite", paired=paired),
			mem=memory,
			cores=1,
			partition=partition
			)

##################################
# function to collect taxonomies #
##################################

def collect():
	command = '''{a}'''.format(
		a = "metaxa2_dc -i " + scratch + "*.level* -o [targets[0]]"
		)
	return str(command)

collect_depends = []
for name in names:
	collect_depends.append(list_depends(name=name, step="collect", paired=paired)[0])

workflow.add_task_gridable(actions=collect(),
	depends=collect_depends,
	targets=list_targets(name=name, step="collect", paired=paired),
	time=calculate_time(name=name, step="collect", paired=paired),
	mem=memory,
	cores=1,
	partition=partition
	)

############################
# merge to MPA-like format #
############################

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

depends = [output + "merged.tsv"]

targets = output + "merged/merged_reorganized.tsv"
if not os.path.isdir(output + "merged/"):
	os.makedirs(output + "merged/")

workflow.add_task("Rscript " + this_folder + "metaxa2_wrangling.R -i [depends[0]] -o [targets[0]]",
	depends = depends,
	targets = targets
	)

####################
# run the workflow #
####################

workflow.go()

#
