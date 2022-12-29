# To run: $ python download_files_grid.py -i for_downloads.txt -o OUTPUT_FOLDER --grid-partition 'shared' --grid-jobs 96 --time 60 --grid-scratch SCRATCH
# In general, download time = 3 GB/hr

from anadama2 import Workflow
import pandas as pd
import glob as glob
import os

workflow = Workflow()

workflow.add_argument("time", desc="The maximum time in minutes allocated to run the command", type=int, default=300)
args = workflow.parse_args()
output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)

scratch = "/" + args.grid_scratch.strip("/") + "/"
if not os.path.isdir(scratch):
	os.makedirs(scratch)

scratch = "/" + args.grid_scratch.strip("/") + "/"
if not os.path.isdir(scratch):
	os.makedirs(scratch)

scratch_tmp = scratch + "tmp/"
if not os.path.isdir(scratch_tmp):
	os.makedirs(scratch_tmp)

df = pd.read_csv(str(args.input), header = None)

files = (df.iloc[:, 0]).tolist()

for file in files:
    paths = glob.glob(output + file + '*.fastq*')
    if len(paths) == 0:
        command = '''{a}'''.format(
            a = "fasterq-dump " + file + " -O " + scratch + " -t " + scratch_tmp + " && if [ -f " + scratch + file + ".fastq ]; then gzip " + scratch + file + ".fastq && cp " + scratch + file + ".fastq.gz " + output + "; else gzip " + scratch + file + "_1.fastq && gzip " + scratch + file + "_2.fastq && cp " + scratch + file + "_1.fastq.gz " + output + " && cp " + scratch + file + "_2.fastq.gz " + output + "; fi;"
            )
        workflow.add_task_gridable(command,
            time=int(args.time),
			mem=4000,
			cores=1,
			partition=args.grid_partition)

workflow.go()
