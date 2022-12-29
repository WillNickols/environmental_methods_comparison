#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 36:00:00
#SBATCH --mem 70000
#SBATCH -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/saltmarsh/saltmarsh_test.out
#SBATCH -e /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/saltmarsh/saltmarsh_test.err
#SBATCH --account=nguyen_lab

mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/
mkdir -p /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/spades_scratch/
mkdir -p /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/spades_testing/saltmarsh/

cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/kneaddata/reERR5004402_paired_1.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/
cp /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/kneaddata/reERR5004402_paired_2.fastq.gz /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/
python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/metaspades/SPAdes-3.15.4-Linux/bin/spades.py --meta -1 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/reERR5004402_paired_1.fastq.gz -2 /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/reERR5004402_paired_2.fastq.gz -o /n/holyscratch01/nguyen_lab/wnickols/spades_testing/saltmarsh/spades_scratch/ -m 70000 -t 16
