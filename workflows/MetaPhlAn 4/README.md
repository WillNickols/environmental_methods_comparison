# MetaPhlAn 4 workflow

This folder contains code for the MetaPhlAn 4 workflow.  The `mpa4_workflow_human.py` script is for when only the `_paired_1.fastq.gz` and `_paired_2.fastq.gz` files are available (this only applied to the human reads).

# Versioning:

- MetaPhlAn4: 4.beta.2 (The most recent version, 4.0.3, includes MAGs generated in this project, so this previous version is used.)
  - Database: mpa_vJan21_CHOCOPhlAnSGB_202103

# Running

Commands of the following form were run:

```
python mpa4_workflow.py -i kneaddata_cleaned \
  -o metaphlan4_output \
  --input-extension fastq.gz \
  --paired paired
```
