#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path

workflow = Workflow(version="0.4", description="mag workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run the command", type=int, default=45000)
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=10000)
workflow.add_argument("phylophlan-time", desc="The maximum time in minutes allocated to run phylophlan", type=int, default=10000)
workflow.add_argument(name="paired", default="paired")
args = workflow.parse_args()
input_extension = args.input_extension

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
paired = args.paired
phylophlan_time = args.phylophlan_time
local_jobs = args.jobs

# list the input fastq files
in_dir = args.input

if paired == "paired":
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)
	files = []
	for path in paths:
		files.append(path)
	names = set(file.split('_paired')[0].split('_unmatched')[0] for file in files)

if paired == "unpaired":
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)
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
		if step == "megahit":
			time = 270 * n_gigabytes
		if step == "metabat":
			time = 30 * n_gigabytes
		if step == "abundance":
			time = 15 * n_gigabytes
	if paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "megahit":
			time = 135 * n_gigabytes
		if step == "metabat":
			time = 15 * n_gigabytes
		if step == "abundance":
			time = 15 * n_gigabytes
	if input_extension == "fastq.gz":
		time = time * 6
	if time > max_time:
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	out_dir = output  + "mags/" + name.split("/")[-1] + "/"
	if paired == "paired":
		if step == "megahit":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
		if step == "metabat":
			depends = out_dir + name.split("/")[-1] + ".contig_depths.txt"
		if step == "abundance":
			depends = out_dir + "bins/" + name.split("/")[-1] + ".bin.unbinned.fa"
	if paired == "unpaired":
		if step == "megahit":
			return [str(name + "." + input_extension)]
		if step == "metabat":
			depends = out_dir + name.split("/")[-1] + ".contig_depths.txt"
		if step == "abundance":
			depends = out_dir + "bins/" + name.split("/")[-1] + ".bin.unbinned.fa"
	return [str(depends)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	out_dir = output + "mags/" + name.split("/")[-1] + "/"
	if step == "megahit":
		targets0 = out_dir + name.split("/")[-1] + ".final.contigs.fa"
		targets1 = out_dir + name.split("/")[-1] + ".contig_depths.txt"
		return [str(targets0), str(targets1)]
	if step == "metabat":
		targets = out_dir + "bins/" + name.split("/")[-1] + ".bin.unbinned.fa"
		return [str(targets)]
	if step == "abundance":
		targets0 = output + "checkm/" + name.split("/")[-1] + ".coverage.tsv"
		targets1 = output + "checkm/" + name.split("/")[-1] + ".taxonomy.tsv"
		targets2 = output + "checkm/" + name.split("/")[-1] + ".mapped_read_num.txt"
		return [str(targets0), str(targets1), str(targets2)]
	return []

###################################################
# function to call MEGAHIT and calculate coverage #
###################################################

def megahit_and_align_and_abundance(name, paired):
	main_dir = scratch + "mags/" + name.split("/")[-1] + "/"
	megahit_dir = main_dir + "megahit/"
	contigs = main_dir + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = scratch + "mags/" + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"

	if input_extension == "fastq.gz":
		if paired == "paired":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h}'''.format(
				a = "megahit -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -r " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -o " + megahit_dir + " -m 0.99 -t " + str(cores) + " --continue --min-contig-len 1500",
				b = "mv " + megahit_dir + "final.contigs.fa " + contigs,
				c = "mkdir -p " + bowtie2_dir,
				d = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
				e = "bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -U " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				f = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				g = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				h = "jgi_summarize_bam_contig_depths --outputDepth [targets[1]] " + bam_sorted,
				)
		if paired == "unpaired":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h}'''.format(
				a = "megahit -r " + name + "." + input_extension + " -o " + megahit_dir + " -m 0.99 -t " + str(cores) + " --continue --min-contig-len 1500",
				b = "mv " + megahit_dir + "final.contigs.fa " + contigs,
				c = "mkdir -p " + bowtie2_dir,
				d = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
				e = "bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				f = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				g = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				h = "jgi_summarize_bam_contig_depths --outputDepth [targets[1]] " + bam_sorted,
				)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="megahit", paired=paired)[0]):
		workflow.add_task_gridable(actions=megahit_and_align_and_abundance(name, paired),
			depends=list_depends(name=name, step="megahit", paired=paired),
			targets=list_targets(name=name, step="megahit", paired=paired),
			time=calculate_time(name=name, step="megahit", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

#############################
# function to run MetaBAT 2 #
#############################

def metabat(name):
	contigs = output + "mags/" + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	depth = output + "mags/" + name.split("/")[-1] + "/" + name.split("/")[-1] + ".contig_depths.txt"
	metabat_tmp = scratch + "mags/" + name.split("/")[-1] + "/bins/" + name.split("/")[-1]
	metabat_out = output + "mags/" + name.split("/")[-1] + "/bins/"
	command = '''{a} && {b} && {c}'''.format(
		a = "if [ ! -s " + contigs + " ]; then touch " + metabat_tmp + ".bin.lowDepth.fa && touch " + metabat_tmp + ".bin.tooShort.fa && touch " + metabat_tmp + ".bin.unbinned.fa; else metabat2 -i " + contigs + " -a " + depth + " -o " + metabat_tmp + ".bin --unbinned -m 1500 -t 3; fi",
		b = "mkdir -p " + metabat_out,
		c = "ls -l " + metabat_tmp + "*.fa | sed 's/.* //g' | grep -v \".bin.unbinned.fa\" | awk '{system (\"cp \" $0 \" " + metabat_out + "\")}'"
		)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="metabat", paired=paired)[0]):
		workflow.add_task_gridable(actions=metabat(name),
			depends=list_depends(name=name, step="metabat", paired=paired),
			targets=list_targets(name=name, step="metabat", paired=paired),
			time=calculate_time(name=name, step="metabat", paired=paired),
			mem=memory,
			cores=1,
			partition=partition
			)

###################################
# function to calculate abundance #
###################################

def abundance(name, paired):
	main_dir = scratch + "mags/" + name.split("/")[-1] + "/"
	megahit_dir = main_dir + "megahit/"
	contigs = main_dir + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = scratch + "mags/" + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"
	bam_index = bowtie2_dir + name.split("/")[-1] + ".sorted.bam.bai"
	bin = output + "mags/" + name.split("/")[-1] + "/bins"

	if paired == "paired":
		if input_extension == "fastq.gz":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r",
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(zcat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				g = "unpaired1=$(echo $(zcat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				h = "unpaired2=$(echo $(zcat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)

	if paired == "unpaired":
		if input_extension == "fastq.gz":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} '''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r",
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "echo $(zcat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				f = "rm -r " + bowtie2_dir
				)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="abundance", paired=paired)[0]):
		workflow.add_task_gridable(actions=abundance(name, paired),
			depends=list_depends(name=name, step="abundance", paired=paired),
			targets=list_targets(name=name, step="abundance", paired=paired),
			time=calculate_time(name=name, step="abundance", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

#######################
# run checkm workflow #
#######################

qa_dir = output + "checkm_qa/"
if not os.path.isdir(qa_dir):
	os.makedirs(qa_dir)

# threads
threads = str(cores)

##################################
# copy bins to to output/checkm/ #
##################################

checkmdir = output+"checkm/"

command = "python checkm_copy_bins.py -a " + output + "mags/" + " -b " + checkmdir + " -n 1000"

workflow.add_task(command, targets=checkmdir, depends=[list_targets(name=name, step="abundance", paired=paired)[0] for name in names])

#########################
# calculate N50 of MAGs #
#########################

if not os.path.exists(output + "n50/mags_n50.tsv"):
	n50 = "python mag_n50_calc.py -i " + checkmdir + " -o " + output + "n50/ -t " + threads
	workflow.add_task(n50, targets=output + "n50/mags_n50.tsv", depends=checkmdir)

##############
# run CheckM #
##############

file=checkmdir
tmp_out = output + "tmp_" + file.split("/")[-1]
qa = qa_dir + file.split("/")[-1] + "checkm_qa.tsv"
if not os.path.exists(qa):
	command = '''{a} && {b}'''.format(
		a = "checkm lineage_wf -x fa " + file + " " + tmp_out + " -f [targets[0]] --tab_table --threads " + str(threads),
		b = "rm -r " + tmp_out
		)
	workflow.add_task(command, depends = [output + "n50/mags_n50.tsv", file], targets=qa)

#######################################
# combine + tidy the CheckM QA tables #
#######################################

if not os.path.exists(output + "merged/checkm_qa_combined.tsv"):
	depends = [output + "n50/mags_n50.tsv"]
	depends.append(qa)

	targets = output + "merged/checkm_qa_combined.tsv"

	if not os.path.isdir(output+"merged/"):
		os.makedirs(output+"merged/")

	workflow.add_task("Rscript checkm_wrangling.R -i [qa_dir] --n50 [n50] -o [targets[0]]",
		depends = depends,
		targets = targets,
		qa_dir = output + "checkm_qa/",
		n50 = output + "n50/mags_n50.tsv"
		)

##################
# run PhyloPhlAn #
##################

if args.grid_options is not None:
	command = "python phylophlan_workflow.py --bins " + output + "mags/" + " --database SGB.Jul20 -o " + output + "phylophlan/ --database-folder /n/holystore01/LABS/huttenhower_lab/Users/wnickols/phylophlan/phylophlan_databases --mem " + str(2 * memory) + " --cores " + str(4 * cores) + " --time " + str(phylophlan_time) + " --input-extension fa --grid-jobs 1 --grid slurm --grid-partition " + partition + " --grid-scratch " + scratch + "phylophlan --local-jobs " + str(local_jobs) + " --grid-options=\"" + " ".join(args.grid_options) + "\" --abundance " + output + "checkm/" + " --qa " + output + "merged/checkm_qa_combined.tsv"
else:
	command = "python phylophlan_workflow.py --bins " + output + "mags/" + " --database SGB.Jul20 -o " + output + "phylophlan/ --database-folder /n/holystore01/LABS/huttenhower_lab/Users/wnickols/phylophlan/phylophlan_databases --mem " + str(2 * memory) + " --cores " + str(4 * cores) + " --time " + str(phylophlan_time) + " --input-extension fa --grid-jobs 1 --grid slurm --grid-partition " + partition + " --grid-scratch " + scratch + "phylophlan --local-jobs " + str(local_jobs) + " --abundance " + output + "checkm/" + " --qa " + output + "merged/checkm_qa_combined.tsv"

workflow.add_task(command, targets = output + "phylophlan/phylophlan_tax_merged.tsv", depends = output + "merged/checkm_qa_combined.tsv")

####################
# run the workflow #
####################

workflow.go()

#
