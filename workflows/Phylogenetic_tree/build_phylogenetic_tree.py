#!/usr/bin/env python

from anadama2 import Workflow
import os
import math

workflow = Workflow(version="0.1", description="Build phylogenetic tree from reference genomes using mashtree")
workflow.add_argument("taxid-list", desc="File with one NCBI TaxID per line")
workflow.add_argument("cores", desc="Number of CPU cores for mashtree", type=int, default=16)
workflow.add_argument("mem", desc="Maximum memory in megabytes", type=int, default=64000)
workflow.add_argument("time", desc="Maximum time in minutes", type=int, default=1440)
workflow.add_argument("sketch-size", desc="Mash sketch size for mashtree", type=int, default=10000)
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

output = os.path.abspath(args.output.rstrip("/")) + "/"
os.makedirs(output, exist_ok=True)

scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/"

taxid_list = os.path.abspath(args.taxid_list)
cores = args.cores
memory = args.mem
partition = args.grid_partition
max_time = args.time
sketch_size = args.sketch_size

genomes_dir = output + "genomes/"
os.makedirs(genomes_dir, exist_ok=True)

output_tree = output + "phylogenetic_tree.nwk"

# Task 1: Download genomes from NCBI FTP (one per species TaxID)
download_done = genomes_dir + "download.done"

download_command = (
    "python {this_folder}download_genomes.py "
    "--taxid-list {taxid_list} "
    "--output-dir {genomes_dir} && "
    "touch {download_done}"
).format(
    this_folder=this_folder,
    taxid_list=taxid_list,
    genomes_dir=genomes_dir,
    download_done=download_done
)

workflow.add_task_gridable(
    actions=download_command,
    depends=[taxid_list],
    targets=[download_done],
    time=max_time,
    mem=memory,
    cores=cores,
    partition=partition
)

# Task 2: Build tree with mashtree
mashtree_command = (
    "mashtree --numcpus {cores} --sketch-size {sketch_size} "
    "{genomes_dir}*.fna > {output_tree}"
).format(
    cores=cores,
    sketch_size=sketch_size,
    genomes_dir=genomes_dir,
    output_tree=output_tree
)

workflow.add_task_gridable(
    actions=mashtree_command,
    depends=[download_done],
    targets=[output_tree],
    time=max_time,
    mem=memory,
    cores=cores,
    partition=partition
)

workflow.go()
