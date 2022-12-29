#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 72:00:00
#SBATCH --mem 184000
#SBATCH -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/human/human_test.out
#SBATCH -e /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/human/human_test.err
#SBATCH --account=nguyen_lab

mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/
mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/spades_scratch/
mkdir -p /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/human/

cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/kneaddata/CSM67UH7_paired_1.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/
cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/kneaddata/CSM67UH7_paired_2.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/
python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/metaspades/SPAdes-3.15.4-Linux/bin/spades.py --meta -1 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/CSM67UH7_paired_1.fastq.gz -2 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/CSM67UH7_paired_2.fastq.gz -o /n/holyscratch01/nguyen_lab/wnickols/spades_testing/human/spades_scratch/ -m 184000 -t 48
