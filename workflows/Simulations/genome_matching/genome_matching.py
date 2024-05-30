#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path
import csv
import pandas as pd

# Run like: 
# conda activate biobakery_assembly
# python genome_matching.py -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/ --cores 1 --abundance-type by_sample --assembler MEGAHIT -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/genome_matching/ --local-jobs 28 --profiler gtdbtk

# Parse arguments
workflow = Workflow(version="0.1", description="MAG and SGB workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("abundance-type", desc='Should abundance estimates come from aligning reads against all MAGs in the run (by_sample) or just the MAGs generated from the sample (by_sample)?', default="by_sample")
workflow.add_argument("assembler", desc='MEGAHIT or metaSPAdes', default="MEGAHIT")
workflow.add_argument("profiler", default="phylophlan")

args = workflow.parse_args()
local_jobs = args.jobs
cores = args.cores

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

# --input should be soil/ ocean/ or animal_gut/
abundance_type = args.abundance_type
try:
	assert abundance_type in ['by_sample', 'by_dataset']
except:
	raise ValueError("--abundance must be by_sample or by_dataset")

input_path = os.path.abspath(args.input.rstrip("/")) + "/"
output_path = os.path.abspath(args.output.rstrip("/")) + "/"
os.makedirs(output_path, exist_ok=True)

sgb_path = input_path + "new_sgbs/sgbs/sgbs/"
genomes_dir = input_path + "simulations/genomes_dir/"

paths = Path(sgb_path).rglob("*.fa")

files = []
for path in paths:
    files.append(str(path))

paths = Path(genomes_dir + "/known/genomes/").rglob("*.fa")

for path in paths:
    files.append(str(path))

paths = Path(genomes_dir + "/new/genomes/").rglob("*.fa")

for path in paths:
    files.append(str(path))

if args.assembler == "MEGAHIT":
	assembly_dir = "assembly"
else:
	assembly_dir = "assembly_metaspades"

mash_dir = output_path + "mash/"
os.makedirs(mash_dir, exist_ok=True)

# List MAGs to run Mash on
with open(mash_dir + 'mags_filepaths.txt', 'w') as f:
    for file in files:
        f.write(file + "\n")

sample_names = []
for potential_dir in os.listdir(input_path + "simulations/" + assembly_dir + "/main/bins/"):
	if not "qc_bins" in potential_dir:
		d = os.path.join(input_path + "simulations/" + assembly_dir + "/main/bins/", potential_dir)
		if os.path.isdir(d):
			paths = Path(d.rstrip("/") + "/bins/").rglob("*.fa")
			files = [str(path) for path in paths]
			files = [file for file in files if not "unbinned" in file and not "tooShort" in file and not "lowDepth" in file and not "final.contigs" in file]

			tmp_df = pd.read_csv(input_path + "simulations/" + assembly_dir + "/main/checkm/qa/checkm_qa_and_n50.tsv", sep="\t")
			bins_to_keep = tmp_df[tmp_df['keep'] == "keep"]['bin_id'].tolist()

			files = [file for file in files if file.split("/")[-1].strip(".fa") in bins_to_keep]

			if len(files) > 0:
				sample_names.append(potential_dir.strip("/"))
				os.makedirs(mash_dir + potential_dir.strip("/"), exist_ok=True)
				with open(mash_dir + potential_dir.strip("/") + '/sample_filepaths.txt', 'w') as f:
					for file in files:
						f.write(file + "\n")

# Mash sketch
sketch = "mash sketch -p " + str(1) + " -l " + mash_dir + "mags_filepaths.txt" + " -o " + mash_dir + "sketches"
workflow.add_task(sketch, depends=mash_dir + "mags_filepaths.txt", targets=mash_dir + "sketches.msh", name="Mash sketch")

for sample_name in sample_names:
	sketch = "mash sketch -p " + str(1) + " -l " + mash_dir + sample_name + '/sample_filepaths.txt' + " -o " + mash_dir + sample_name + "/sketches"
	workflow.add_task(sketch, depends=mash_dir + sample_name + '/sample_filepaths.txt', targets=mash_dir + sample_name + "/sketches.msh", name="Mash sketch for " + sample_name)

	dist = "mash dist -p " + str(1) + " -t " + mash_dir + sample_name + "/sketches.msh " + mash_dir + "sketches.msh > " + mash_dir + sample_name + "/mash_dist_out.tsv"
	workflow.add_task(dist, depends=[mash_dir + sample_name + "/sketches.msh", mash_dir + "sketches.msh"], targets=mash_dir + sample_name + "/mash_dist_out.tsv", name="Mash dist")

matching = "Rscript " + this_folder + "parse_and_merge.R -i " + input_path + "simulations/" + assembly_dir + "/main/" + " --abundance_dir " + input_path + "simulations/" + assembly_dir + "/main/abundance_" + abundance_type + "/" + " --mash " + mash_dir + " -o " + output_path + "genome_matching.tsv" + " --profiler " + args.profiler
depends_list = [mash_dir + sample_name + "/mash_dist_out.tsv" for sample_name in sample_names]

if args.profiler == "phylophlan":
	depends_list.append(input_path + "simulations/" + assembly_dir + "/main/phylophlan/phylophlan_relab.tsv")
else:
	depends_list.append(input_path + "simulations/" + "gtdbtk/relab/gtdbtk_relab.tsv")
depends_list.append(input_path + "simulations/" + assembly_dir + "/main/checkm/qa/checkm_qa_and_n50.tsv")
workflow.add_task(matching, depends=depends_list, targets=output_path + "genome_matching.tsv")

####################
# run the workflow #
####################

workflow.go()
