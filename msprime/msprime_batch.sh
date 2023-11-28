#!/bin/bash -l
#SBATCH -J aqu_pol
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/msprime/logs/%a.out
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-100
#SBATCH --ntasks 1
#SBATCH --mem=2G

export PATH="/projappl/project_2001443/msprime-env/bin:$PATH"

mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/sim_$SLURM_ARRAY_TASK_ID
cd  /scratch/project_2001443/barriers_introgr_formica/msprime/sim_$SLURM_ARRAY_TASK_ID

python3.6 /scratch/project_2001443/barriers_introgr_formica/msprime/aqu_pol_msprime.py
