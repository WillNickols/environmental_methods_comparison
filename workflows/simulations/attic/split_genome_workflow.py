# To run: $ python split_genome_workflow.py -i $INPUT_DIRECTORY -o $OUTPUT_DIRECTORY split_length retain_length

from anadama2 import Workflow
import os
import glob
import pandas as pd

workflow = Workflow()

workflow.add_argument("split-length", desc="split length", default = 100)
workflow.add_argument("retain-length", desc="retain length", default = 1000)
workflow.add_argument("metadata", desc="metadata file")
args = workflow.parse_args()

output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
    os.makedirs(output)

paths = glob.glob("/" + args.input.strip("/") + "/*")

metadata = pd.read_csv(args.metadata, sep='\t')

def splitGenome(file):
	command = '''{a} && {b} && {c}'''.format(
        a = "gunzip " + file,
        b = "python splitGenome.py " + file.split(".gz")[0] + " " + output + (file.split("/")[-1]).split("_genomic")[0] + "_" + args.split_length + "_" + args.retain_length + ".fa " + args.split_length + " " + args.retain_length,
        c = "gzip " + file.split(".gz")[0],
		)
	return str(command)

for i in range(metadata.shape[0]):
    strain = metadata['strain_ID'][i]
    if not os.path.isfile(output + strain + "_" + args.split_length + "_" + args.retain_length + ".fa"):
        workflow.add_task(actions=splitGenome(metadata['filepath'][i]), depends=[metadata['filepath'][i]], targets=[output + strain + "_" + args.split_length + "_" + args.retain_length + ".fa"])

    metadata['filepath'][i] = output + strain + "_" + args.split_length + "_" + args.retain_length + ".fa"

workflow.go()

metadata.to_csv(args.metadata, sep = "\t", index = False)
