#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 72:00:00
#SBATCH --mem 184000
#SBATCH -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/wild_gut/wild_gut_test.out
#SBATCH -e /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/wild_gut/wild_gut_test.err
#SBATCH --account=nguyen_lab

mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/wild_gut/
mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/wild_gut/spades_scratch/
mkdir -p /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/wild_gut/

cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/reERR5005154_paired_1.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/wild_gut/
cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/reERR5005154_paired_2.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/wild_gut/
python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/metaspades/SPAdes-3.15.4-Linux/bin/spades.py --meta -1 /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/reERR5005154_paired_1.fastq.gz -2 /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/reERR5005154_paired_2.fastq.gz -o /n/holyscratch01/nguyen_lab/wnickols/spades_testing/wild_gut/spades_scratch/ -m 184000 -t 48
