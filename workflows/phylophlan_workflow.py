#!/usr/bin/env python

from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
from pathlib import Path
import glob
import os
import math
import shutil

workflow = Workflow(version="0.4", description="CheckM abundance workflow")
workflow.add_argument("bins", desc="Directory with all bins from mags workflow", default="mags/")
workflow.add_argument("database", desc="Database name", default="SGB.Jul20")
workflow.add_argument("database-folder", desc="Folder with the phylophlan database")
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=20000)
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=6)
workflow.add_argument(name="input-extension", desc="The bin file extensions", default="fastq")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("abundance", desc="checkm_tax abundance folder")
workflow.add_argument("qa", desc="checkm_qa_combined file")

args = workflow.parse_args()
input_extension = args.input_extension

# bins
bins = "/" + args.bins.strip("/") + "/"

paths = Path("/" + bins.strip("/") + "/").rglob('*/bins/*.' + input_extension)

files = []
not_allowed = ["unbinned", "tooShort", "lowDepth", "final.contigs"]
for path in paths:
	if not any(x in str(path) for x in not_allowed):
		files.append(path.as_posix())

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

# database
database = args.database
database_folder = args.database_folder

#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(fastq, step):
	if step == "copy_bins":
		time = 5
	if step == "phylophlan":
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(fastq, step):
	if step == "copy_bins":
		depends = fastq
	if step == "phylophlan":
		depends = [str(scratch + "bins/" + fastq.split("/")[-1]) for fastq in files]
		return depends
	return [str(depends)]

############################
# function to list targets #
############################

def list_targets(fastq, step):
	full_name = fastq.split("/")[-1]
	name = fastq.split("." + input_extension)[0].split("/")[-1]
	if step == "copy_bins":
		targets = scratch + "bins/" + full_name
	if step == "phylophlan":
		targets = output + "combined.tsv"
		return [str(targets)]
	return [str(targets)]

###########################
# run abundance profiling #
###########################

def copy_bins(fastq):
	phylophlan_dir = scratch + "bins/"
	command = '''{a} && {b}'''.format(
		a = "mkdir -p " + phylophlan_dir,
		b = "cp " + fastq + " " + phylophlan_dir
		)
	return str(command)

for file in files:
	if not os.path.isfile(list_targets(fastq=file, step="copy_bins")[0]):
		workflow.add_task(actions=copy_bins(file),
			depends=list_depends(fastq=file, step="copy_bins"),
			targets=list_targets(fastq=file, step="copy_bins"))

def phylophlan(fastq):
	name = fastq.split("." + input_extension)[0].split("/")[-1]
	phylophlan_dir = scratch + "bins/"
	output_name = scratch + "combined"

	command = '''{a}'''.format(
		a = "phylophlan_metagenomic -i " + phylophlan_dir + " -n 1 --add_ggb --add_fgb -d " + database + " -o " + output_name + " --nproc " + str(cores) + " --verbose --database_folder " + database_folder,
		)
	return str(command)

if not os.path.isfile(list_targets(fastq=file, step="phylophlan")[0]):
	workflow.add_task_gridable(actions=phylophlan(file),
		depends=[str(scratch + "bins/" + fastq.split("/")[-1]) for fastq in files],
		targets=output + "combined.tsv",
		time=max_time,
		mem=memory,
		cores=cores,
		partition=partition
		)

wrangling_dep = [args.abundance, output+"combined.tsv", args.qa]
wrangling_targets = [output+"phylophlan_tax_merged.tsv"]

workflow.add_task("Rscript phylophlan_tax_wrangling.R -i [depends[0]] --tax [depends[1]] --qa [depends[2]] -o [targets[0]]",
	depends = wrangling_dep,
	targets = wrangling_targets
	)

################
# run workflow #
################

workflow.go()

#
