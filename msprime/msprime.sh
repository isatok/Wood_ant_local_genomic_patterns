#---------------------------------------
# FINNISH AQU, FINNISH POL
#---------------------------------------

#!/bin/bash -l
#SBATCH -J aquFI_polFI
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/msprime/logs/aquFI_polFI_%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/msprime/logs/aquFI_polFI_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1
#SBATCH --ntasks 1
#SBATCH --mem=2G

#change later on #SBATCH --array=1-100 if replicates are needed

export PATH="/projappl/project_2001443/msprime-env/bin:$PATH"

mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/aquFI_polFI_sim_$SLURM_ARRAY_TASK_ID
cd  /scratch/project_2001443/barriers_introgr_formica/msprime/aquFI_polFI_sim_$SLURM_ARRAY_TASK_ID

python /scratch/project_2001443/barriers_introgr_formica/msprime/aquFI_polFI_msprime.py

###END.

#---------------------------------------
# SWISS AQU, WEST-SWISS POL
#---------------------------------------


#!/bin/bash -l
#SBATCH -J aquSWI_polWSWI
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/msprime/logs/aquSWI_polWSWI_%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/msprime/logs/aquSWI_polWSWI_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1
#SBATCH --ntasks 1
#SBATCH --mem=2G

#change later on #SBATCH --array=1-100 if replicates are needed

export PATH="/projappl/project_2001443/msprime-env/bin:$PATH"

mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/aquSWI_polWSWI_sim_$SLURM_ARRAY_TASK_ID
cd  /scratch/project_2001443/barriers_introgr_formica/msprime/aquSWI_polWSWI_sim_$SLURM_ARRAY_TASK_ID

python /scratch/project_2001443/barriers_introgr_formica/msprime/aquSWI_polWSWI_msprime.py

###END.
