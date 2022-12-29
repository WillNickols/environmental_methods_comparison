from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import shutil
from pathlib import Path

# python split_inputs_workflow.py -i INPUT_FOLDER -o OUTPUT_FOLDER --input-extension fastq.gz

workflow = Workflow(version="0.4", description="split interleaved files for kneaddata")
workflow.add_argument(name="input-extension", desc="the input file extension", default="fastq.gz")
args = workflow.parse_args()

paths = Path(os.path.abspath(args.input.rstrip("/")) + "/").rglob('*.' + str(args.input_extension))

files = []
for path in paths:
	files.append(path.as_posix())

output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

if str(args.input_extension) == "fastq.gz":
    out_files = [output + file.split("/")[-1].split("." + args.input_extension)[0] for file in files]

for f,o in zip(files, out_files):
	if not os.path.isfile(o + "_1.fastq.gz") or not os.path.isfile(o + "_2.fastq.gz"):
    	workflow.add_task("bash reformat.sh in=" + f + " out=" + o + "_1.fastq.gz out2=" + o + "_2.fastq.gz ", depends = [f], targets = [o + "_1.fastq.gz", o + "_2.fastq.gz"])
workflow.go()
