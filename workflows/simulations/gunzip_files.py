from anadama2 import Workflow

workflow = Workflow()

args = workflow.parse_args()
in_files = workflow.get_input_files(input_folder = args.input, extension=".fna.gz")
out_files = [file.strip(".fna.gz") + ".fa" for file in in_files]

for in_file, out_file in zip(in_files, out_files):
    workflow.add_task(
        "gunzip [depends[0]] -f && mv " + in_file.strip(".gz") + " " + out_file,
        depends=in_file,
        targets=out_file)

workflow.go()
