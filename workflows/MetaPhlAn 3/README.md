# MetaPhlAn 3 workflow

This folder contains code for the MetaPhlAn 3 workflow.  The `mpa3_workflow_human.py` script is for when only the `_paired_1.fastq.gz` and `_paired_2.fastq.gz` files are available (this only applied to the human reads).

# Versioning:

- MetaPhlAn3: 3.0.14
  - Database: mpa_v30_CHOCOPhlAn_201901

# Running

Commands of the following form were run:

```
python mpa3_workflow.py -i kneaddata_cleaned \
  -o metaphlan3_output \
  --input-extension fastq.gz \
  --paired paired
```
