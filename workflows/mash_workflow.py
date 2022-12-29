#!/usr/bin/env python

from anadama2 import Workflow
unused_options = ["input", "grid_jobs", "grid", "grid_partition", "grid_benchmark", "grid_options", "grid_environment", "grid_scratch"]

workflow = Workflow(version="0.1", description="Mash workflow", remove_options=unused_options)
workflow.add_argument("checkm-qa", desc="CheckM QA table", type=str, required=True)
workflow.add_argument("bins", desc="Directory containing MAGs", type=str, required=True)
workflow.add_argument("threads", desc="The number of threads to run per task", type=int, default=1)
args = workflow.parse_args()

checkm = args.checkm_qa
bins = "/" + args.bins.strip("/") + "/"
threads = str(args.threads)
out_dir = "/" + args.output.strip("/") + "/"

mash_dir = out_dir + "mash/"

############################
# list MAGs to run Mash on #
############################

list_inputs = "python mash_list_inputs.py --checkm " + checkm + " --bins " + bins + " --mash " + mash_dir + " --threads " + threads

workflow.add_task(list_inputs, targets=mash_dir + "mags_filepaths.txt")

###############
# mash sketch #
###############

sketch = "mash sketch -p " + threads + " -l " + mash_dir + "mags_filepaths.txt" + " -o " + mash_dir + "sketches"

workflow.add_task(sketch, depends=mash_dir + "mags_filepaths.txt", targets=mash_dir + "sketches.msh")

##############
# mash paste #
##############

paste = "mash paste " + mash_dir + "references " + mash_dir + "sketches.msh"

workflow.add_task(paste, depends=mash_dir + "sketches.msh", targets=mash_dir + "references.msh")

#############
# mash dist #
#############

dist = "mash dist -p " + threads + " -t " + mash_dir + "references.msh " + mash_dir + "sketches.msh > " + mash_dir + "mash_dist_out.tsv"

workflow.add_task(dist, depends=mash_dir + "references.msh", targets=mash_dir + "mash_dist_out.tsv")

#################
# identify SGBS #
#################

# groups MAGs into SGBs using (1) Mash and (2) fastANI

cluster = "Rscript mash_clusters.R --mash " + mash_dir + "mash_dist_out.tsv --checkm " + checkm + " --out_dir " + out_dir + " --threads " + threads + " --mag_dir " + bins + "qc_bins/"

workflow.add_task(cluster, depends=mash_dir + "mash_dist_out.tsv", targets=[out_dir + "fastANI/SGB_list.txt", out_dir + "sgbs/SGB_info.tsv"])

####################
# run the workflow #
####################

workflow.go()
