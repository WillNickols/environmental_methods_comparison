# To run: $ python generateUnknownSequences.py --local-jobs 4 -o $OUTPUT_FOLDER --start-date YYYY/MM/DD --end-date YYYY/MM/DD --n --split-length --retain-length

from anadama2 import Workflow
import os

workflow = Workflow(version = "0.4")

workflow.add_argument("start-date", desc="start date to download sequences", default = "2022/05/11")
workflow.add_argument("end-date", desc="end date to download sequences", default = "2022/06/11")
workflow.add_argument("n", desc="number of sequences to download", default = 10)
workflow.add_argument("split-length", desc="split length", default = 100)
workflow.add_argument("retain-length", desc="retain length", default = 1000)
workflow.add_argument("tax-dump-folder", desc="NCBI tax dump")
args = workflow.parse_args()

output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
    os.makedirs(output)

tax_dump_folder = args.tax_dump_folder

workflow.add_task(
    "python downloadGenomesForShredding_ncbi_update.py " + args.start_date + " " + args.end_date + " " + str(args.n) + " " + output + " " + tax_dump_folder,
    targets=output + "shred_metadata.tsv")

workflow.add_task(
    "python split_genome_workflow.py -i " + output + "genomes/ -o " + output + "genomes/ --split-length " + str(args.split_length) + " --retain-length " + str(args.retain_length) + " --local-jobs " + str(args.jobs) + " --metadata " + output + "shred_metadata.tsv",
    depends=output + "shred_metadata.tsv")

workflow.go()
