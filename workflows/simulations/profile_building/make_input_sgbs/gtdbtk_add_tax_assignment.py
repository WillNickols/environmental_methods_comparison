import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--bac", help="gtdbtk bac120 table", type=str, required=True)
parser.add_argument("--arc", help="gtdbtk ar53 table", type=str, required=True)
parser.add_argument("--output", help="output file", type=str, required=True)
args = parser.parse_args()

bac = pd.read_csv(args.bac, sep='\t')
arc = pd.read_csv(args.arc, sep='\t')

table = pd.concat([bac, arc], axis=0)
table = table[['user_genome', 'classification']]
table.rename(columns={"user_genome": "mag", "classification": "taxon"}, inplace=True)
table.loc[table['taxon'].isin(["Unclassified", "Unclassified Bacteria", "Unclassified Archaea"]), 'taxon'] = "UNKNOWN"
table['taxon'] = table['taxon'].str.replace(';','|')
table['taxon'] = table['taxon'].str.replace('d__','k__')

table.to_csv(args.output, sep='\t', index=False)
