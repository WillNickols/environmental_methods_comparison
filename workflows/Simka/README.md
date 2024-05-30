# Simka workflow

This folder contains code for the Simka comparison workflow.

# Versioning:

- Simka: 1.5.1

# Installation

Simka was installed according to the instructions on [its GitHub](https://github.com/GATB/simka).

# Running

Commands of the following form were run:

```
python simka_workflow.py \
  --input kneaddata_cleaned \
  --output simka_output \
  --paired paired \
  --simka-path simka/build/bin/simka
```

```
python simka_workflow_all.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simka_all/all_kneaddata.txt \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simka_all/output/ \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/simka_all/ \
  --grid-partition 'shared' --grid-jobs 100 --cores 40 --mem 184000 --local-jobs 8 \
  --grid-options="--account=nguyen_lab"
```