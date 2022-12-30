# mOTUs3 workflow

This folder contains code for the mOTUs3 workflow.

# Versioning:

- mOTUs3: 3.0.1
  - Database: db_mOTU_v3.0.1

# Installation

mOTUs3 was installed from Conda according to [the GitHub instructions](https://github.com/motu-tool/mOTUs).

# Running

Commands of the following form were run:

```
python mOTUs3_workflow.py -i kneaddata_cleaned \
  -o mOTUs3_output \
  --input-extension fastq.gz \
  --paired paired
```
