# Pre-processing scripts

This folder contains code to download read files and clean them.  A previous version of Kneaddata had a known bug with SRR and ERR read headers, so the rehead script makes the reads compatible with all Kneaddata versions.  Kneaddata was run with `--serial` and `--run-trf`, and the contaminant database was built from the hg19 human genome except for on the dog gut and cat gut samples.  The dog gut samples used a Kneaddata database built from the canis lupus familiaris genome (GCF_014441545), and the cat samples used a Kneaddata database built from the felis catus genome (GCA_000181335). 

# Versioning

- fasterq-dump: 2.11.0
- Kneaddata: 0.12.0

# Running

Commands of the following form were run:

```
python download_files_grid.py -i files_to_download.tsv \
  -o raw_inputs/

python rehead_workflow.py -i raw_inputs/ \
  -o reheaded/ \
  --input-extension fastq.gz

python kneaddata_workflow.py -i reheaded/ \
  -o kneaddata_cleaned/ \
  --input-extension fastq.gz \
  --paired paired
```