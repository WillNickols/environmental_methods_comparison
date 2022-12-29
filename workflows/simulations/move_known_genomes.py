import pandas as pd
import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--general-reference-metadata", help="general_reference_metadata.tsv")
parser.add_argument("--reference-profile", help="reference_profile.tsv")
parser.add_argument("--src", help="directory to copy from")
parser.add_argument("--dest", help="directory to copy to")
args = parser.parse_args()

gen_ref = pd.read_csv(args.general_reference_metadata, sep = '\t')
df = pd.read_csv(args.reference_profile, sep = '\t')

column_names = ["genome_ID", "NCBI_ID", "filepath", "novelty_category",  "strain_ID"]
output_df = pd.DataFrame(columns = column_names)

for tax_id in df['NCBI_ID']:
    new_rows = gen_ref.loc[gen_ref['NCBI_ID'] == tax_id]
    if len(new_rows.index) > 0:
        for file in new_rows['filepath'].tolist():
            shutil.copyfile(file, file.replace(args.src, args.dest))
        new_rows['filepath'] = new_rows['filepath'].str.replace(args.src, args.dest)
        output_df = pd.concat([output_df,new_rows], ignore_index = True)

output_df.to_csv(args.dest.rsplit("genomes/", 1)[0] + "/reference_metadata.tsv", index = False, sep = "\t")
