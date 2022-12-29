#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 120:00:00
#SBATCH --mem 184000
#SBATCH -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/tara/tara_test.out
#SBATCH -e /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/tara/tara_test.err
#SBATCH --account=nguyen_lab

mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/
mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/spades_scratch/
mkdir -p /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/tara/

cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/kneaddata/reERR3589580_paired_1.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/
cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/kneaddata/reERR3589580_paired_2.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/
python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/metaspades/SPAdes-3.15.4-Linux/bin/spades.py --meta -1 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/reERR3589580_paired_1.fastq.gz -2 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/reERR3589580_paired_2.fastq.gz -o /n/holyscratch01/nguyen_lab/wnickols/spades_testing/tara/spades_scratch/ -m 184000 -t 16
