#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path

# Parse arguments
workflow = Workflow(version="0.1", description="Spades workflow using contigs as inputs")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run any individual command except metaSPADES", type=int, default=180000)
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run any individual command", type=int, default=10000)
workflow.add_argument("time-per-gb", desc="Time in hr per gigabyte for metaSPADES run")
workflow.add_argument("mem-per-gb", desc="Memory in GB per gigabyte for metaSPADES run")
workflow.add_argument("paired", desc="Whether the inputs are \"paired\", or something else", default="paired")
workflow.add_argument("pair-identifier", desc="extension to identify the first file in a pair like _R1; \"kneaddata_default\" looks for _paired_1, _paired_2 by default", default="kneaddata_default")
workflow.add_argument("skip-placement", desc="Whether to stop after checkm steps", action="store_true")
workflow.add_argument("remove-intermediate-output", desc="Remove intermediate files", action="store_true")
workflow.add_argument("min-contig-length", desc='MetaBAT -m parameter', default=1500)
workflow.add_argument("metabat-options", desc='MetaBAT options as a text string with quotes', default="")
workflow.add_argument("checkm-coverage-options", desc='checkm coverage options as a text string with quotes', default="")
workflow.add_argument("checkm-predict-options", desc='checkm2 predict options as a text string with quotes', default="")
workflow.add_argument("abundance-type", desc='Should abundance estimates come from aligning reads against all MAGs in the run (by_dataset) or just the MAGs generated from the sample (by_sample)?', default="by_sample")
workflow.add_argument("completeness", desc='completeness threshold for retaining bins', default=50)
workflow.add_argument("contamination", desc='contamination threshold for retaining bins', default=10)
workflow.add_argument("checkm-predict-options", desc='checkm2 predict options as a text string with quotes', default="")
workflow.add_argument("phylophlan-metagenomic-options", desc='PhyloPhlAn metagenomic options as a text string with quotes', default="")
workflow.add_argument("gc-length-stats", desc='calculate GC and length stats for each bin', action="store_true")
workflow.add_argument("mash-sketch-options", desc='Mash sketch options as a text string with quotes', default="")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

assembly_tasks_folder = this_folder + "assembly_tasks/"

try:
	os.environ["CHECKM_DATA_PATH"]
except:
	raise ValueError("CHECKM_DATA_PATH not provided or set")

if not args.skip_placement:
	try:
		database_folder = os.environ["PHYLOPHLAN_PATH"]
		database_files = os.listdir(os.environ["PHYLOPHLAN_PATH"])
		database = [file for file in database_files if file.count('.') == 1][0]
	except:
		raise ValueError("PHYLOPHLAN_PATH not provided or set")

# Check valid input extension
input_extension = args.input_extension
try:
	assert input_extension in ['fastq', 'fastq.gz', 'fq', 'fq.gz']
except:
	raise ValueError("--input-extension must be fastq, fastq.gz, fq, or fq.gz")

paired = args.paired
try:
	assert paired in ['paired']
except:
	raise ValueError("--paired must be paired")
pair_identifier = args.pair_identifier
pair_identifier_2 = pair_identifier.replace("1", "2")

abundance_type = args.abundance_type
try:
	assert abundance_type in ['by_sample', 'by_dataset']
except:
	raise ValueError("--abundance must be by_sample or by_dataset")

# list the input fastq files
in_dir = args.input

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time
local_jobs = args.jobs

#######################################
# function to calculate tool runtimes #
#######################################

time_per_gb = float(args.time_per_gb)
mem_per_gb = float(args.mem_per_gb)

def calculate_time(name, step, no_max = False):
	time = 0
	if step == "abundance_dataset":
		if pair_identifier == "kneaddata_default":
			time = 50 * sum([math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0)) for name in names])
		else:
			time = 50 * sum([math.ceil(os.path.getsize(name + pair_identifier + "." + input_extension) / (1024 * 1024 * 1024.0)) for name in names])
		time = time + 30
	else:
		if pair_identifier == "kneaddata_default":
			n_gigabytes = os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0)
		else:
			n_gigabytes = os.path.getsize(name + pair_identifier + "." + input_extension) / (1024 * 1024 * 1024.0)
		if step == "spades":
			time = time_per_gb * 60 * n_gigabytes + 3
		elif step == "align":
			time = 600 * n_gigabytes + 3
		elif step == "metabat":
			time = 600 * n_gigabytes + 3
		elif step == "abundance":
			time = 400 * n_gigabytes + 3
		elif step == "checkm":
			time = 300 * n_gigabytes + 3
		elif step == "phylophlan":
			time = 300 * n_gigabytes + 20
	if 'gz' in input_extension:
		time = time * 6
	if not no_max:
		if time > max_time:
			time = max_time
	return math.ceil(time)

# Get filepath without paired ending
paths = glob.glob(os.path.abspath(in_dir.rstrip("/")) + "/" + '*.' + input_extension)
if pair_identifier == "kneaddata_default":
	names_tmp = set(file.split('_paired_1')[0].split('_paired_2')[0].split('_unmatched_1')[0].split('_unmatched_2')[0] for file in paths)
else:
	names_tmp = set(file.split(pair_identifier)[0].split(pair_identifier_2)[0] for file in paths)

if pair_identifier == "kneaddata_default":
	names = [name for name in names_tmp if calculate_time(name, "spades", True) < args.time and math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024.0) * mem_per_gb) < 184000]
else:
	names = [name for name in names_tmp if calculate_time(name, "spades", True) < args.time and math.ceil(os.path.getsize(name + pair_identifier + "." + input_extension) / (1024 * 1024.0) * mem_per_gb) < 184000]

# output
original_output = os.path.abspath(args.output.rstrip("/")) + "/"
os.makedirs(original_output, exist_ok=True)

output = original_output + "main/"
os.makedirs(output, exist_ok=True)

# scratch directory
if args.grid_scratch:
	scratch = os.path.abspath(args.grid_scratch.rstrip("/")) + "/main/"
else:
	scratch = output + "scratch/"
os.makedirs(scratch, exist_ok=True)

assembly_dir = output + "assembly/"
os.makedirs(assembly_dir, exist_ok=True)

spades_scratch = scratch + "spades/"
os.makedirs(spades_scratch, exist_ok=True)

contigs_dir = assembly_dir + "main/"
os.makedirs(contigs_dir, exist_ok=True)

names = [name for name in names if os.path.exists(str(contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"))]

print(set(names_tmp).difference(set(names)))

if len(names) == 0:
	raise ValueError("No input files")

mags_scratch = scratch + "bins/"
os.makedirs(mags_scratch, exist_ok=True)

bowtie2_tmp = scratch + "bowtie2_tmp/"
os.makedirs(bowtie2_tmp, exist_ok=True)

bowtie2_global_dir = scratch + "bowtie2_global_dir/"
os.makedirs(bowtie2_global_dir, exist_ok=True)

depths_dir = assembly_dir + "contig_depths/"
os.makedirs(depths_dir, exist_ok=True)

bins_dir = output + "bins/"
os.makedirs(bins_dir, exist_ok=True)

abundance_dir = output + "abundance_" + abundance_type + "/"
os.makedirs(abundance_dir, exist_ok=True)

checkm_dir = output + "checkm/"
os.makedirs(checkm_dir, exist_ok=True)

checkm_scratch = scratch + "checkm/"
os.makedirs(checkm_scratch, exist_ok=True)

checkm_qa_scratch = checkm_scratch + "qa_unmerged/"
os.makedirs(checkm_qa_scratch, exist_ok=True)

checkm_bins_dir_scratch = checkm_scratch + "bins/"
os.makedirs(checkm_bins_dir_scratch, exist_ok=True)

checkm_n50_dir = checkm_dir + "n50/"
os.makedirs(checkm_n50_dir, exist_ok=True)

qa_unmerged_dir = checkm_dir + "qa_unmerged/"
os.makedirs(qa_unmerged_dir, exist_ok=True)

qa_dir = checkm_dir + "qa/"
os.makedirs(qa_dir, exist_ok=True)

if not args.skip_placement:
	phylophlan_scratch = scratch + "phylophlan/"
	os.makedirs(phylophlan_scratch, exist_ok=True)

	phylophlan_dir = output + "phylophlan/"
	os.makedirs(phylophlan_dir, exist_ok=True)

	phylophlan_unmerged_dir = phylophlan_dir + "unmerged/"
	os.makedirs(phylophlan_unmerged_dir, exist_ok=True)

	phylophlan_unmerged_scratch = phylophlan_scratch + "unmerged/"
	os.makedirs(phylophlan_unmerged_scratch, exist_ok=True)

	sgb_dir = output + "sgbs/"
	os.makedirs(sgb_dir, exist_ok=True)

	mash_dir = sgb_dir + "mash/"
	os.makedirs(mash_dir, exist_ok=True)

	sgbs_scratch = scratch + "sgbs/"
	os.makedirs(sgbs_scratch, exist_ok=True)

	sgbs = sgb_dir + "sgbs/"

	map_out = sgbs_scratch + "abundances/"
	os.makedirs(map_out, exist_ok=True)

#################################
# function to list dependencies #
#################################

def list_depends(name, step):
	if step == "spades":
		if pair_identifier == "kneaddata_default":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
		else:
			return [str(name + pair_identifier + "." + input_extension), str(name + pair_identifier_2 + "." + input_extension)]
	elif step == "align":
		depends_list = [str(contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa")]
		return depends_list
	elif step == "metabat":
		depends_list = [str(depths_dir + name.split("/")[-1] + ".contig_depths.txt")]
		return depends_list
	elif step == "abundance_sample":
		if pair_identifier == "kneaddata_default":
			depends_list = [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension)]
		else:
			depends_list = [str(name + pair_identifier + "." + input_extension), str(name + pair_identifier_2 + "." + input_extension)]
		depends_list.append(str(bins_dir + name.split("/")[-1] + ".done"))
		return depends_list
	elif step == "copy_bins":
		return [str(bins_dir + name.split("/")[-1] + ".done")]


############################
# function to list targets #
############################

def list_targets(name, step):
	if step == "spades":
		return [str(contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa")]
	elif step == "align":
		return [str(depths_dir + name.split("/")[-1] + ".contig_depths.txt")]
	elif step == "metabat":
		return [str(bins_dir + name.split("/")[-1] + ".done")]
	elif step == "abundance":
		targets0 = abundance_dir + name.split("/")[-1] + ".coverage.tsv"
		targets1 = abundance_dir + name.split("/")[-1] + ".abundance.tsv"
		targets2 = abundance_dir + name.split("/")[-1] + ".mapped_read_num.txt"
		return [str(targets0), str(targets1), str(targets2)]
	elif step == "copy_bins":
		return[str(checkm_bins_dir_scratch + name.split("/")[-1] + "/bins/" + name.split("/")[-1] + ".done")]

###############################
# function to call metaSPADES #
###############################

def spades(name):
	if pair_identifier == "kneaddata_default":
		n_gigabytes = os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0)
	else:
		n_gigabytes = os.path.getsize(name + pair_identifier + "." + input_extension) / (1024 * 1024 * 1024.0)
	if pair_identifier == "kneaddata_default":
		command = '''{a} && {b} && {c}'''.format(
			a = os.environ["SPADES_BIN"] + "spades.py --meta -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -o " + spades_scratch + name.split("/")[-1] + "/ -m " + str(math.ceil(mem_per_gb * n_gigabytes)) + " -t " + str(cores) + " " + args.spades_options,
			b = "mkdir -p " + contigs_dir + name.split("/")[-1],
			c = "mv " + spades_scratch + name.split("/")[-1] + "/contigs.fasta " + contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
			)
	else:
		command = '''{a} && {b} && {c}'''.format(
			a = os.environ["SPADES_BIN"] + "spades.py --meta -1 " + name + pair_identifier + "." + input_extension + " -2 " + name + pair_identifier_2 + "." + input_extension + " -o " + spades_scratch + name.split("/")[-1] + "/ -m " + str(math.ceil(mem_per_gb * n_gigabytes)) + " -t " + str(cores) + " " + args.spades_options,
			b = "mkdir -p " + contigs_dir + name.split("/")[-1],
			c = "mv " + spades_scratch + name.split("/")[-1] + "/contigs.fasta " + contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
			)
	return str(command)

for name in names:
	if not os.path.exists(list_targets(name=name, step="spades")[0]):
		if pair_identifier == "kneaddata_default":
			n_gigabytes = os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024.0)
		else:
			n_gigabytes = os.path.getsize(name + pair_identifier + "." + input_extension) / (1024 * 1024.0)
		workflow.add_task_gridable(actions=spades(name),
			depends=list_depends(name=name, step="spades"),
			targets=list_targets(name=name, step="spades"),
			time=calculate_time(name=name, step="spades"),
			mem=math.ceil(mem_per_gb * n_gigabytes),
			cores=cores,
			partition=partition,
			name="Assemble contigs with spades on " + name.split("/")[-1]
			)

###############################################
# function to align reads and calculate depth #
###############################################

def align(name):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"

	if pair_identifier == "kneaddata_default":
		command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			a = "mkdir -p " + bowtie2_dir,
			b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
			c = "bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
			d = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
			e = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
			f = "jgi_summarize_bam_contig_depths --outputDepth " + depths_dir + name.split("/")[-1] + ".contig_depths.txt " + bam_sorted,
			)
	else:
		command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			a = "mkdir -p " + bowtie2_dir,
			b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
			c = "bowtie2 -x " + index + " -1 " + name + pair_identifier + "." + input_extension + " -2 " + name + pair_identifier_2 + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
			d = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
			e = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
			f = "jgi_summarize_bam_contig_depths --outputDepth " + depths_dir + name.split("/")[-1] + ".contig_depths.txt " + bam_sorted,
			)

	return str(command)

for name in names:
	workflow.add_task_gridable(actions=align(name),
		depends=list_depends(name=name, step="align"),
		targets=list_targets(name=name, step="align"),
		time=calculate_time(name=name, step="align"),
		mem=memory,
		cores=cores,
		partition=partition,
		name="Bowtie2 alignment of " + name.split("/")[-1]
		)

#############################
# function to run MetaBAT 2 #
#############################

def metabat(name):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	depth = depths_dir + name.split("/")[-1] + ".contig_depths.txt"
	metabat_tmp = mags_scratch + name.split("/")[-1] + "/bins/" + name.split("/")[-1]
	metabat_out = bins_dir + name.split("/")[-1] + "/bins/"
	command = '''{a} && {b} && {c} && {d} && {e}'''.format(
		a = "if [ ! -s " + contigs + " ]; then mkdir -p " + metabat_tmp + " && touch " + metabat_tmp + ".bin.lowDepth.fa && touch " + metabat_tmp + ".bin.tooShort.fa && touch " + metabat_tmp + ".bin.unbinned.fa; else metabat2 -i " + contigs + " -a " + depth + " -o " + metabat_tmp + ".bin --unbinned -m 1500 -t " + str(cores) + " " + args.metabat_options + "; fi",
		b = "if [ ! \"$(ls -A " + mags_scratch + name.split("/")[-1] + "/bins/" + ")\" ]; then > " + contigs + " && touch " + metabat_tmp + ".bin.lowDepth.fa && touch " + metabat_tmp + ".bin.tooShort.fa && touch " + metabat_tmp + ".bin.unbinned.fa; fi", # if contigs are too short and no bins are created, delete contigs from the original file and create the empty bin files
		c = "mkdir -p " + metabat_out,
		d = "cp " + metabat_tmp + "*.fa " + metabat_out,
		e = "touch [targets[0]]"
		)
	return str(command)

for name in names:
	workflow.add_task_gridable(actions=metabat(name),
		depends=list_depends(name=name, step="metabat"),
		targets=list_targets(name=name, step="metabat"),
		time=calculate_time(name=name, step="metabat"),
		mem=memory,
		cores=cores,
		partition=partition,
		name="Bin with MetaBAT on " + name.split("/")[-1]
		)

###################################
# function to calculate abundance #
###################################

def abundance_sample(name):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"
	bam_index = bowtie2_dir + name.split("/")[-1] + ".sorted.bam.bai"
	bin = bins_dir + name.split("/")[-1] + "/bins"

	if input_extension in ["fastq.gz", "fq.gz"]:
		if pair_identifier == "kneaddata_default":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(zcat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
		else:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(zcat " + name + pair_identifier + "." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(zcat " + name + pair_identifier_2 + "." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
	else:
		if pair_identifier == "kneaddata_default":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(cat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(cat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
		else:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(cat " + name + pair_identifier + "." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(cat " + name + pair_identifier_2 + "." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
	return str(command)

def rebuild_bowtie2_db():
	command = '''{a} && {b} && {c} && {d}'''.format(
		a = "rm -r " + bowtie2_global_dir,
		b = "mkdir -p " + bowtie2_global_dir,
		c = "python " + assembly_tasks_folder + "by_dataset_db.py --mag_dir " + bins_dir + " --out_dir " + bowtie2_global_dir,
		d = "touch " + abundance_dir + "built.done"
		)
	return str(command)

def abundance_dataset(name):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_global_dir + "mag_db"
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"
	bam_index = bowtie2_dir + name.split("/")[-1] + ".sorted.bam.bai"
	bin = bowtie2_global_dir + "tmp/"

	if input_extension in ["fastq.gz", "fq.gz"]:
		if pair_identifier == "kneaddata_default":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(zcat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2)) &>> [targets[2]]",
				)
		else:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + pair_identifier + "." + input_extension + " -2 " + name + pair_identifier_2 + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(zcat " + name + pair_identifier + "." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(zcat " + name + pair_identifier_2 + "." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2)) &>> [targets[2]]",
				)
	else:
		if pair_identifier == "kneaddata_default":
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(cat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(cat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2)) &>> [targets[2]]",
				)
		else:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + pair_identifier + "." + input_extension + " -2 " + name + pair_identifier_2 + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "python " + assembly_tasks_folder + "checkm.py coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "python " + assembly_tasks_folder + "checkm.py profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(cat " + name + pair_identifier + "." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(cat " + name + pair_identifier_2 + "." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2)) &>> [targets[2]]",
				)
	return str(command)

if abundance_type == "by_sample":
	for name in names:
		workflow.add_task_gridable(actions=abundance_sample(name),
			depends=list_depends(name=name, step="abundance_sample"),
			targets=list_targets(name=name, step="abundance"),
			time=calculate_time(name=name, step="abundance"),
			mem=memory,
			cores=cores,
			partition=partition,
			name="Calculate by-sample abundance for " + name.split("/")[-1]
			)
else:
	depends_list = [list_depends(name, "abundance_sample") for name in names]
	workflow.add_task_gridable(actions=rebuild_bowtie2_db(),
		depends=[item for sublist in depends_list for item in sublist],
		targets=abundance_dir + "built.done",
		time=calculate_time(name="", step="abundance_dataset"),
		mem=memory,
		cores=1,
		partition=partition,
		name="Build bowtie2 database for by-dataset abundance"
		)
	for name in names:
		workflow.add_task_gridable(actions=abundance_dataset(name),
			depends=abundance_dir + "built.done",
			targets=list_targets(name=name, step="abundance"),
			time=calculate_time(name=name, step="abundance"),
			mem=memory,
			cores=cores,
			partition=partition,
			name="Calculate by-dataset abundance for " + name.split("/")[-1]
			)

################
# Calculate GC #
################

if args.gc_length_stats and not os.path.isdir(qa_dir + "GC_content.tsv"):
	command = "python " + assembly_tasks_folder + "calculateGC.py --in-dir " + bins_dir + " --out-file " + qa_dir + "GC_content.tsv" + " --threads " + str(local_jobs)
	workflow.add_task(command, targets=qa_dir + "GC_content.tsv", depends=[item for name in names for item in list_targets(name=name, step="abundance")], name="Calculate GC statistics")

######################################
# run checkm and phylophlan workflow #
######################################

def copy_bins(name):
	metabat_out = bins_dir + name.split("/")[-1] + "/bins/"
	checkm_bin_name = checkm_bins_dir_scratch + name.split("/")[-1] + "/bins/"
	os.makedirs(checkm_bin_name, exist_ok=True)
	command = '''{a} && {b}'''.format(
		a = "if ls " + metabat_out + "*.bin.[0-9]*.fa; then cp " + metabat_out + "*.bin.[0-9]*.fa " + checkm_bin_name + "; fi",
		b = "touch [targets[0]]"
		)
	return str(command)

for name in names:
	workflow.add_task(actions=copy_bins(name),
		depends=list_depends(name=name, step="copy_bins"),
		targets=list_targets(name=name, step="copy_bins"),
		name="Copy bins for checkm for " + name.split("/")[-1]
		)

# Calculate n50 of MAGs
n50 = "python " + assembly_tasks_folder + "mag_n50_calc.py -i " + checkm_bins_dir_scratch + " -o " + checkm_n50_dir + " -t " + str(cores)
workflow.add_task(n50, targets=checkm_n50_dir + "mags_n50.tsv", depends=[list_targets(name=name, step="copy_bins")[0] for name in names], name="Calculate n50")

# Run CheckM
for name in names:
	checkm_bin_name = checkm_bins_dir_scratch + name.split("/")[-1] + "/bins/"
	metabat_out = bins_dir + name.split("/")[-1] + "/bins/"
	os.makedirs(checkm_qa_scratch + name.split("/")[-1] + "/", exist_ok=True)
	os.makedirs(qa_unmerged_dir + name.split("/")[-1] + "/", exist_ok=True)
	if args.grid_jobs > 0:
		command = "if ls " + metabat_out + "*.bin.[0-9]*.fa; then " + this_folder + "checkm2/bin/checkm2 predict -x fa --input " + checkm_bin_name + " --output-directory " + checkm_qa_scratch + name.split("/")[-1] + "/" + " --threads " + str(args.cores) + " " + args.checkm_predict_options + "; " + \
			"else echo -e \"Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes\" > " + qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv; fi"
	else:
		command = "if ls " + metabat_out + "*.bin.[0-9]*.fa; then " + this_folder + "checkm2/bin/checkm2 predict -x fa --input " + checkm_bin_name + " --output-directory " + qa_unmerged_dir + name.split("/")[-1] + "/" + " --threads " + str(args.cores) + " " + args.checkm_predict_options + "; " + \
			"else echo -e \"Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tAdditional_Notes\" > " + qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv; fi"
	workflow.add_task_gridable(actions=command,
	depends=list_targets(name=name, step="copy_bins"),
	targets=qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv",
	time=calculate_time(name, "checkm"),
	mem=memory,
	cores=args.cores,
	partition=partition,
	name="Assess bin quality with checkm2 for " + name.split("/")[-1])

command = ""
depends_list = []
for count, name in enumerate(names):
	if count == 0:
		command = "cat " + qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv > " + qa_dir + "quality_report.tsv"
		depends_list.append(qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv")
	else:
		command = command + " && tail -n +2 " + qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv >> " + qa_dir + "quality_report.tsv"
		depends_list.append(qa_unmerged_dir + name.split("/")[-1] + "/quality_report.tsv")
workflow.add_task(command,
depends = depends_list,
targets = qa_dir + "quality_report.tsv",
name="Merge checkm2 quality report"
)

workflow.add_task("python " + assembly_tasks_folder + "checkm_wrangling.py --checkm-qa " + qa_dir + "quality_report.tsv" + " --n50 " + checkm_n50_dir + "mags_n50.tsv" + " --out_file " + qa_dir + "checkm_qa_and_n50.tsv" + " --completeness " + str(args.completeness) + " --contamination " + str(args.contamination),
	depends = [checkm_n50_dir + "mags_n50.tsv", qa_dir + "quality_report.tsv"],
	targets = qa_dir + "checkm_qa_and_n50.tsv",
	name="Merge bin statistics"
	)

if not args.skip_placement:
	##################
	# run PhyloPhlAn #
	##################
	for name in names:
		os.makedirs(phylophlan_unmerged_scratch + name.split("/")[-1] + "/", exist_ok=True)
		os.makedirs(phylophlan_unmerged_dir + name.split("/")[-1] + "/", exist_ok=True)
		metabat_out = bins_dir + name.split("/")[-1] + "/bins/"
		checkm_bin_name = checkm_bins_dir_scratch + name.split("/")[-1] + "/bins/"
		if args.grid_jobs > 0:
			command = "if ls " + metabat_out + "*.bin.[0-9]*.fa; then phylophlan_metagenomic -i " + checkm_bin_name + " -n 1 --add_ggb --add_fgb -d " + database + " -o " + phylophlan_unmerged_scratch + name.split("/")[-1] + "/" + "phylophlan_out --nproc " + str(cores) + " --verbose -e fa --database_folder " + database_folder + " " + args.phylophlan_metagenomic_options + "; " + \
				"else echo -e \"line1\nline2\nline3\n#mag\tsgb\tggb\tfgb\tref\" > " + phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv; fi"
		else:
			command = "if ls " + metabat_out + "*.bin.[0-9]*.fa; then phylophlan_metagenomic -i " + checkm_bin_name + " -n 1 --add_ggb --add_fgb -d " + database + " -o " + phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out --nproc " + str(cores) + " --verbose -e fa --database_folder " + database_folder + " " + args.phylophlan_metagenomic_options + "; " + \
				"else echo -e \"line1\nline2\nline3\n#mag\tsgb\tggb\tfgb\tref\" > " + phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv; fi"
		workflow.add_task_gridable(actions=command,
		depends=list_targets(name=name, step="copy_bins"),
		targets=phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv",
		time=calculate_time(name, "phylophlan"),
		mem=args.mem,
		cores=args.cores,
		partition=partition,
		name="Run PhyloPhlAn metagenomic on " + name.split("/")[-1])

	command = ""
	depends_list = []
	for count, name in enumerate(names):
		if count == 0:
			command = "cat " + phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv > " + phylophlan_dir + "phylophlan_out.tsv"
			depends_list.append(phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv")
		else:
			command = command + " && tail -n +5 " + phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv >> " + phylophlan_dir + "phylophlan_out.tsv"
			depends_list.append(phylophlan_unmerged_dir + name.split("/")[-1] + "/" + "phylophlan_out.tsv")
	workflow.add_task(command,
	depends = depends_list,
	targets = phylophlan_dir + "phylophlan_out.tsv",
	name="Merge PhyloPhlAn quality report"
	)

	workflow.add_task("python " + assembly_tasks_folder + "phylophlan_add_tax_assignment.py --table " + phylophlan_dir + "phylophlan_out.tsv" + " --output " + phylophlan_dir + "phylophlan_relab.tsv",
	depends=phylophlan_dir + "phylophlan_out.tsv",
	targets=phylophlan_dir + "phylophlan_relab.tsv",
	name="Relabel PhyloPhlAn output")

	########
	# SGBs #
	########

	############################
	# list MAGs to run Mash on #
	############################

	list_inputs = "python " + assembly_tasks_folder + "mash_list_inputs.py --checkm " + qa_dir + "checkm_qa_and_n50.tsv --phylophlan " + phylophlan_dir + "phylophlan_relab.tsv" + " --bins " + bins_dir + " --mash " + mash_dir + " --threads " + str(local_jobs)
	workflow.add_task(list_inputs, depends=[phylophlan_dir + "phylophlan_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv"], targets=mash_dir + "mags_filepaths.txt", name="List Mash inputs")

	###############
	# mash sketch #
	###############

	sketch = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "sketches.msh" + "; else mash sketch -p " + str(local_jobs) + " -l " + mash_dir + "mags_filepaths.txt" + " -o " + mash_dir + "sketches " + args.mash_sketch_options + "; fi"
	workflow.add_task(sketch, depends=mash_dir + "mags_filepaths.txt", targets=mash_dir + "sketches.msh", name="Mash sketch")

	##############
	# mash paste #
	##############

	paste = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "references.msh" + "; else mash paste " + mash_dir + "references " + mash_dir + "sketches.msh; fi"
	workflow.add_task(paste, depends=mash_dir + "sketches.msh", targets=mash_dir + "references.msh", name="Mash paste")

	#############
	# mash dist #
	#############

	dist = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "mash_dist_out.tsv" + "; else mash dist -p " + str(local_jobs) + " -t " + mash_dir + "references.msh " + mash_dir + "sketches.msh > " + mash_dir + "mash_dist_out.tsv; fi"
	workflow.add_task(dist, depends=mash_dir + "references.msh", targets=mash_dir + "mash_dist_out.tsv", name="Mash dist")

	#################
	# identify SGBS #
	#################

	depends_list = [item for name in names for item in list_targets(name=name, step="abundance")]
	depends_list.extend([phylophlan_dir + "phylophlan_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv", sgb_dir + "sgbs/SGB_info.tsv"])

	# groups MAGs into SGBs using (1) Mash and (2) fastANI
	cluster = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then mkdir -p " + sgb_dir + "fastANI/" + " && touch " + sgb_dir + "fastANI/SGB_list.txt && mkdir -p " + sgb_dir + "sgbs/" + " && echo -e \"cluster\tgenome\tcluster_members\tn_genomes\tcompleteness\tcontamination\tstrain_heterogeneity\tn50\tquality\tkeep\tcluster_name\tsgb\" > " + sgb_dir + "sgbs/SGB_info.tsv; else Rscript " + assembly_tasks_folder + "mash_clusters.R --mash " + mash_dir + "mash_dist_out.tsv --checkm " + qa_dir + "checkm_qa_and_n50.tsv" + " --phylo " + phylophlan_dir + "phylophlan_relab.tsv --out_dir " + sgb_dir + " --threads " + str(local_jobs) + " --mag_dir " + bins_dir + "qc_bins/; fi"
	workflow.add_task(cluster, depends=[mash_dir + "mash_dist_out.tsv", phylophlan_dir + "phylophlan_relab.tsv"], targets=[sgb_dir + "sgbs/SGB_info.tsv"], name="Create SGBs")

	merge = "Rscript " + assembly_tasks_folder + "merge_tax_and_abundance.R" + " -i " + abundance_dir + " --tax " + phylophlan_dir + "phylophlan_relab.tsv" + " --qa " + qa_dir + "checkm_qa_and_n50.tsv --sgbs " + sgb_dir + "sgbs/SGB_info.tsv" + " -o " + original_output + "final_profile_" + abundance_type + ".tsv"
	workflow.add_task(merge, depends=depends_list, targets = original_output + "final_profile_" + abundance_type + ".tsv", name="Merge abundance and taxonomy for final output")

if args.remove_intermediate_output:
	rm_command = "rm -r " + scratch + " && touch " + output + "remove_scratch.done"
	depends_list = [qa_dir + "checkm_qa_and_n50.tsv"]
	if not args.skip_placement:
		depends_list.append(original_output + "final_profile_" + abundance_type + ".tsv")
	workflow.add_task(rm_command, depends=depends_list, targets=output + "remove_scratch.done", name="Remove scratch files")

####################
# run the workflow #
####################

workflow.go()
