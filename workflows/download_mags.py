# To run: $ python download_files.py -i for_downloads.txt -o OUTPUT_FOLDER
# In general, download time = 3 GB/hr

from anadama2 import Workflow
import pandas as pd
import glob as glob
import os

workflow = Workflow()

args = workflow.parse_args()
output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)

df = pd.read_csv(str(args.input), header = None)

files = (df.iloc[:, 0]).tolist()

for file in files:
    paths = glob.glob(output + file + '*.fasta')
    if len(paths) == 0:
        command = '''{a} && {b} && {c}'''.format(
            a = "wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/" + file[0:2] + "/" + file[2:4] + "/" + file[4:6] + "/" + file + "/" + file + ".1.fsa_nt.gz -P " + output,
			b = "gunzip " + output + file + ".1.fsa_nt.gz",
			c = "mv " + output + file + ".1.fsa_nt " + output + file + ".fasta"
            )
        workflow.add_task(command)

workflow.go()
