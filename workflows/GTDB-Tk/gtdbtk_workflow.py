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
workflow.add_argument("mags-sgb-dir", desc='mags and sgb directory')
workflow.add_argument("mash-sketch-options", desc='Mash sketch options as a text string with quotes', default="")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

output = os.path.abspath(args.output.rstrip("/")) + "/"
bins_dir = os.path.abspath(args.input.rstrip("/")) + "/"
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

qa_dir = args.mags_sgb_dir + "checkm/qa/"
abundance_dir = args.mags_sgb_dir + "abundance_" + abundance_type + "/"

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

########
# SGBs #
########

############################
# list MAGs to run Mash on #
############################

if not os.path.isfile(mash_dir + "mags_filepaths.txt"):
	list_inputs = "python " + this_folder + "mash_list_inputs.py --checkm " + qa_dir + "checkm_qa_and_n50.tsv --gtdbtk-table " + gtdbtk_relab + "gtdbtk_relab.tsv" + " --bins " + bins_dir + " --mash " + mash_dir + " --threads " + str(local_jobs)
	workflow.add_task(list_inputs, depends=[gtdbtk_relab + "gtdbtk_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv"], targets=mash_dir + "mags_filepaths.txt")

###############
# mash sketch #
###############

if not os.path.isfile(mash_dir + "sketches.msh"):
	sketch = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "sketches.msh" + "; else mash sketch -p " + str(local_jobs) + " -l " + mash_dir + "mags_filepaths.txt" + " -o " + mash_dir + "sketches " + args.mash_sketch_options + "; fi"
	workflow.add_task(sketch, depends=mash_dir + "mags_filepaths.txt", targets=mash_dir + "sketches.msh")

##############
# mash paste #
##############

if not os.path.isfile(mash_dir + "references.msh"):
	paste = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "references.msh" + "; else mash paste " + mash_dir + "references " + mash_dir + "sketches.msh; fi"
	workflow.add_task(paste, depends=mash_dir + "sketches.msh", targets=mash_dir + "references.msh")

#############
# mash dist #
#############

if not os.path.isfile(mash_dir + "mash_dist_out.tsv"):
	dist = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "mash_dist_out.tsv" + "; else mash dist -p " + str(local_jobs) + " -t " + mash_dir + "references.msh " + mash_dir + "sketches.msh > " + mash_dir + "mash_dist_out.tsv; fi"
	workflow.add_task(dist, depends=mash_dir + "references.msh", targets=mash_dir + "mash_dist_out.tsv")

#################
# identify SGBS #
#################

depends_list = [gtdbtk_relab + "gtdbtk_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv", sgb_dir + "sgbs/SGB_info.tsv"]

# groups MAGs into SGBs using (1) Mash and (2) fastANI
if not os.path.isfile(sgb_dir + "fastANI/SGB_list.txt") or not os.path.isfile(sgb_dir + "sgbs/SGB_info.tsv"):
	cluster = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then mkdir -p " + sgb_dir + "fastANI/" + " && touch " + sgb_dir + "fastANI/SGB_list.txt && mkdir -p " + sgb_dir + "sgbs/" + " && echo -e \"cluster\tgenome\tcluster_members\tn_genomes\tcompleteness\tcontamination\tstrain_heterogeneity\tn50\tquality\tkeep\tcluster_name\tsgb\" > " + sgb_dir + "sgbs/SGB_info.tsv; else Rscript " + this_folder + "mash_clusters.R --mash " + mash_dir + "mash_dist_out.tsv --checkm " + qa_dir + "checkm_qa_and_n50.tsv" + " --gtdbtk " + gtdbtk_relab + "gtdbtk_relab.tsv --out_dir " + sgb_dir + " --threads " + str(local_jobs) + " --mag_dir " + bins_dir + "qc_bins/; fi"
	workflow.add_task(cluster, depends=[mash_dir + "mash_dist_out.tsv", gtdbtk_relab + "gtdbtk_relab.tsv"], targets=[sgb_dir + "fastANI/SGB_list.txt", sgb_dir + "sgbs/SGB_info.tsv"])

if not os.path.isfile(output + "final_profile.tsv"):
	merge = "Rscript " + this_folder + "merge_tax_and_abundance.R" + " -i " + abundance_dir + " --tax " + gtdbtk_relab + "gtdbtk_relab.tsv" + " --qa " + qa_dir + "checkm_qa_and_n50.tsv --sgbs " + sgb_dir + "sgbs/SGB_info.tsv" + " -o " + output + "final_profile_" + abundance_type + ".tsv"
	workflow.add_task(merge, depends=depends_list, targets = output + "final_profile_" + abundance_type + ".tsv")

workflow.go()
