#!/usr/bin/env python

from anadama2 import Workflow
from pathlib import Path
import math
import os
import subprocess
import gzip

workflow = Workflow(version="0.4", description="Kmer distance workflow")

workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in gigabytes allocated to run the command", type=int, default=16000)
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("simka-path", desc="folder with build/bin/simka")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)
input = "/" + args.input.strip("/") + "/"
scratch = "/" + args.grid_scratch.strip("/") + "/"
if not os.path.isdir(scratch):
	os.makedirs(scratch)
simka_path = "/" + args.simka_path.strip("/") + "/"
cores = args.cores
memory = args.mem
partition = args.grid_partition
max_time = args.time

paths = Path(input).rglob('*.fastq.gz')

names = [str(path).split("/")[-1].split("_paired")[0].split("_unmatched")[0].split(".fastq.gz")[0] for path in paths]

time = 15
for i, name in enumerate(names):
	with open(output + 'file_list.txt', 'a') as file_list:
		if os.path.isfile(input + name + '_paired_1.fastq.gz'):
			line_to_write = 'ID' + str(i) + ': ' + input + name + '_paired_1.fastq.gz'
			with gzip.open(input + name + '_unmatched_1.fastq.gz', 'rb') as f:
				if f.readline() != b'':
					line_to_write = line_to_write + ', ' + input + name + '_unmatched_1.fastq.gz'
			with gzip.open(input + name + '_unmatched_2.fastq.gz', 'rb') as f:
				if f.readline() != b'':
					line_to_write = line_to_write + ', ' + input + name + '_unmatched_2.fastq.gz; '
				else:
					line_to_write = line_to_write + '; '
			line_to_write = line_to_write + ' ' + input + name + '_paired_2.fastq.gz\n'
			file_list.write(line_to_write)
			time += 3 * math.ceil(os.path.getsize(input + name + '_paired_1.fastq.gz') / (1024 * 1024 * 1024.0))
		else:
			file_list.write('ID' + str(i) + ': ' + input + name + '.fastq.gz\n')
			time += math.ceil(os.path.getsize(input + name + '.fastq.gz') / (1024 * 1024 * 1024.0))

subprocess.run("cp " + output + 'file_list.txt ' + scratch, shell=True)

# run simka

command = '''{a} && {b} && {c}'''.format(
	a = simka_path + "build/bin/simka -in " + scratch + "file_list.txt" + " -out " + scratch + " -out-tmp " + scratch + "tmp/ -nb-cores " + str(cores) + " -max-memory " + str(memory),
	b = "gunzip " + output + "mat_abundance_braycurtis.csv.gz",
	c = "gunzip " + output + "mat_presenceAbsence_jaccard.csv.gz"
	)

if not os.path.isfile(output + "mat_presenceAbsence_jaccard.csv") or not os.path.isfile(output + "mat_abundance_braycurtis.csv"):
	workflow.add_task_gridable(actions=command,
		depends=scratch + "file_list.txt",
		targets=[output + "mat_presenceAbsence_jaccard.csv", output + "mat_abundance_braycurtis.csv"],
		time=time,
		mem=memory,
		cores=cores,
		partition=partition
		)

workflow.go()
