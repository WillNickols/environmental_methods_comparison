#!/usr/bin/env python

from anadama2 import Workflow
from pathlib import Path
import math
import os

workflow = Workflow(version="0.4", description="GTDBTK abundance workflow")

workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in gigabytes allocated to run the command", type=int, default=16000)
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("abundance-type", desc='by_sample or by_dataset', default="by_sample")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

output = os.path.abspath(args.output.rstrip("/")) + "/"
input = os.path.abspath(args.input.rstrip("/")) + "/"
scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"
cores = args.cores
memory = args.mem
partition = args.grid_partition
max_time = args.time
local_jobs = args.jobs

abundance_type = args.abundance_type
try:
	assert abundance_type in ['by_sample', 'by_dataset']
except:
	raise ValueError("--abundance must be by_sample or by_dataset")

bins_dir = input + "bins/"
paths = Path(bins_dir).rglob('*.fa')

bins_scratch = scratch + "bins_in/"
if not os.path.isdir(bins_scratch):
	os.makedirs(bins_scratch)

main_scratch = scratch + "main/"
if not os.path.isdir(main_scratch):
	os.makedirs(main_scratch)

gtdbtk_out = output + "raw/"
if not os.path.isdir(gtdbtk_out):
	os.makedirs(gtdbtk_out)

gtdbtk_relab = output + "relab/"
if not os.path.isdir(gtdbtk_relab):
	os.makedirs(gtdbtk_relab)

mash_dir = output + "mash/"
if not os.path.isdir(mash_dir):
	os.makedirs(mash_dir)

sgb_dir = output + "sgbs/"
if not os.path.isdir(sgb_dir):
	os.makedirs(sgb_dir)

qa_dir = input + "checkm/qa/"
abundance_dir = input + "abundance_" + abundance_type + "/"

##############
# run GTDBTK #
##############

command = '''{a} && {b} && {c} && {d}'''.format(
	a = "find " + bins_dir + " ! '(' -name '*unbinned*' -o -name '*lowDepth*' -o -name '*tooShort*' -o -name '*final.contigs*' -o -name '*.done' -o -path '*qc_bins*' ')' -type f -exec cp -t " + bins_scratch + " {} +",
	b = "gtdbtk classify_wf --genome_dir " + bins_scratch + " --out_dir " + scratch + "raw/" + " -x .fa --cpus " + str(cores)  + " --scratch_dir " + main_scratch,
	c = "if [ ! -f " + scratch + "raw/" + "gtdbtk.ar53.summary.tsv" + " ]; then echo -e \"user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\" > " + scratch + "raw/" + "gtdbtk.ar53.summary.tsv; fi",
	d = "if [ ! -f " + scratch + "raw/" + "gtdbtk.bac120.summary.tsv" + " ]; then echo -e \"user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\" > " + scratch + "raw/" + "gtdbtk.bac120.summary.tsv; fi"
	)

if not os.path.isfile(gtdbtk_out + "gtdbtk.bac120.summary.tsv") or not os.path.isfile(gtdbtk_out + "gtdbtk.ar53.summary.tsv"):
	workflow.add_task_gridable(actions=command,
		depends=bins_dir,
		targets=[gtdbtk_out + "gtdbtk.bac120.summary.tsv", gtdbtk_out + "gtdbtk.ar53.summary.tsv"],
		time=min(200 + 10 * len(list(paths)), max_time),
		mem=memory,
		cores=cores,
		partition=partition
		)

if not os.path.isfile(gtdbtk_relab + "gtdbtk_relab.tsv"):
	workflow.add_task("python " + this_folder + "gtdbtk_add_tax_assignment.py --bac " + gtdbtk_out + "gtdbtk.bac120.summary.tsv" + " --arc " + gtdbtk_out + "gtdbtk.ar53.summary.tsv --output " + gtdbtk_relab + "gtdbtk_relab.tsv",
		depends=[gtdbtk_out + "gtdbtk.bac120.summary.tsv", gtdbtk_out + "gtdbtk.ar53.summary.tsv"],
		targets=gtdbtk_relab + "gtdbtk_relab.tsv")

depends_list = [gtdbtk_relab + "gtdbtk_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv"]

if not os.path.isfile(output + "final_profile.tsv"):
	merge = "Rscript " + this_folder + "merge_tax_and_abundance.R" + " -i " + abundance_dir + " --tax " + gtdbtk_relab + "gtdbtk_relab.tsv" + " --qa " + qa_dir + "checkm_qa_and_n50.tsv -o " + output + "final_profile_" + abundance_type + ".tsv"
	workflow.add_task(merge, depends=depends_list, targets = output + "final_profile_" + abundance_type + ".tsv")

workflow.go()
