import pandas as pd
from Bio import SeqIO
import os

assembly_summary = pd.read_csv("library/archaea/assembly_summary.txt", sep="\t", skiprows=2, header=None)
assembly_summary = pd.concat([assembly_summary, pd.read_csv("library/bacteria/assembly_summary.txt", sep="\t", skiprows=2, header=None)], ignore_index=True)
assembly_summary = pd.concat([assembly_summary, pd.read_csv("library/viral/assembly_summary.txt", sep="\t", skiprows=2, header=None)], ignore_index=True)

seqid2taxid = pd.read_csv("seqid2taxid.map", sep="    ", header=None)

genomes_to_remove = set(assembly_summary[(assembly_summary.iloc[:,14] > '2022/03/25')].iloc[:,0].tolist())

for domain in ['archaea', 'bacteria', 'viral']:
	assembly_summary = pd.read_csv("library/" + domain + "/assembly_summary.txt", sep="\t", skiprows=2, header=None)
	assembly_summary_filtered = pd.read_csv("library/" + domain + "/assembly_summary_filtered.txt", sep="\t", skiprows=2, header=None)
	for genome in os.listdir('library/' + domain):
		if genome.split('_')[0] + '_' + genome.split('_')[1] in genomes_to_remove:
			fasta_sequences = SeqIO.parse(open('library/' + domain + '/' + genome),'fasta')
			fasta_names = [fasta.id for fasta in fasta_sequences]

			# Remove sequences corresponding to genome from seqid2taxid.map
			seqid2taxid = seqid2taxid[~seqid2taxid.iloc[:,0].isin(fasta_names)]

			# Remove genome from library
			print('Removing ' + 'library/' + domain + '/' + genome)
			os.remove('library/' + domain + '/' + genome)

			assembly_summary = assembly_summary[~assembly_summary.iloc[:,0].isin([genome.split('_')[0] + '_' + genome.split('_')[1]])]
			assembly_summary_filtered = assembly_summary_filtered[~assembly_summary_filtered.iloc[:,0].isin([genome.split('_')[0] + '_' + genome.split('_')[1]])]

	with open("library/" + domain + "/assembly_summary.txt", 'r') as f1:
		with open("library/" + domain + "/assembly_summary_tmp.txt", 'w') as f2:
			for line in f1.readlines()[0:2]:
				f2.write(line)

	assembly_summary.to_csv("library/" + domain + "/assembly_summary_tmp.txt", mode='a', index=False, header=False, sep="\t")
	assembly_summary_filtered.to_csv("library/" + domain + "/assembly_summary_filtered_tmp.txt", index=False, header=False, sep="\t")

with open("seqid2taxid_tmp.map", 'w') as f:
	for seqid, taxid in zip(seqid2taxid.iloc[:,0].tolist(), seqid2taxid.iloc[:,1].tolist()):
		f.write(str(seqid) + "    " + str(taxid) + "\n")
