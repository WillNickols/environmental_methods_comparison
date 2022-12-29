# To run: $ python gzip_files.gz -i $FOLDER --local-jobs 4 -o $OUTPUT

from anadama2 import Workflow

workflow = Workflow()

args = workflow.parse_args()
in_files = workflow.get_input_files(input_folder = args.input, extension=".fastq")
out_files = [file+".gz"  for file in in_files]

workflow.add_task_group(
    "gzip [depends[0]]",
    depends=in_files,
    targets=out_files)

workflow.go()
