# To run: $ python download_files.py -i for_downloads.txt -o OUTPUT_FOLDER
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

df = pd.read_csv(str(args.input), header = None)

files = (df.iloc[:, 0]).tolist()

for file in files:
    paths = glob.glob(output + file + '*.fastq*')
    if len(paths) == 0:
        command = '''{a}'''.format(
            a = "fasterq-dump " + file + " -O " + output + " && if [ -f " + output + file + ".fastq ]; then gzip " + output + file + ".fastq; else gzip " + output + file + "_1.fastq && gzip " + output + file + "_2.fastq; fi;"
            )
        workflow.add_task(command)

workflow.go()
