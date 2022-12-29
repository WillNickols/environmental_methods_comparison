from anadama2 import Workflow
from anadama2.tracked import TrackedDirectory
import glob
import os
import shutil
from pathlib import Path

# python rehead_workflow.py -i INPUT_FOLDER -o OUTPUT_FOLDER --input-extension fastq.gz --local-jobs 16

workflow = Workflow(version="0.1", description="fix fastq sequence headers")
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
    out_files = [output + "re" + (file.split("/")[-1]).split(".gz")[0] for file in files]
elif str(args.input_extension) == "fastq":
    out_files = [output + "re" + file.split("/")[-1] for file in files]

for f,o in zip(files, out_files):
    workflow.add_task("python rehead.py " + f + " " + o + " && gzip " + o, depends = [f], targets = [o + ".gz"])
workflow.go()
