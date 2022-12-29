import gzip
import sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later

# Use: python rehead.py input_filename output_filename

reheaded = open(str(sys.argv[2]), "w")

name = (str(sys.argv[2])).split(".fastq")[0]
pair_id = ""
if name[-2:] == "_1":
    pair_id = "/1"
elif name[-2:] == "_2":
    pair_id = "/2"

if ".gz" in str(sys.argv[1]):
    with gzip.open(str(sys.argv[1]), "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            title = title.rsplit(" ")[0] + pair_id
            reheaded.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
else:
    with open(str(sys.argv[1]), "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            title = title.rsplit(" ")[0] + pair_id
            reheaded.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

reheaded.close()
