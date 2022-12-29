import numpy as np
import sys

# Use: python splitGenome.py input output split_length retain_length

input = str(sys.argv[1])
output = str(sys.argv[2])
split_length = int(sys.argv[3])
retain_length = int(sys.argv[4])

with open(input) as f:
    fasta = "".join([line.strip() for line in f if ">" not in line])

indices_start = np.arange(0, len(fasta), split_length + retain_length).tolist()
indices_end = [index + retain_length for index in indices_start]
indices_end.pop()
indices_end.append(None)
parts = [fasta[i:j] for i,j in zip(indices_start, indices_end)]

with open(output, 'w') as f:
    for i in range(len(parts)):
        f.write(">Simulated_read_" + str(i + 1) + "\n")
        f.write(parts[i] + "\n")
