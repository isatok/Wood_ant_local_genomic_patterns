#### This script is to compute local Pi, Dxy and Fst estimates with PIXY for simulated aqu-pol data ####

mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/pixy
mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/pixy/groupfiles
mkdir /scratch/project_2001443/barriers_introgr_formica/msprime/pixy/logs

cd /scratch/project_2001443/barriers_introgr_formica/msprime/pixy


###
### 1. Prepare a group file for PIXY ------
###

#for simulated data; see file group_pixy_sim.tab
tsk_0	aqu
tsk_1	aqu
tsk_2	aqu
tsk_3	aqu
tsk_4	aqu
tsk_5	aqu
tsk_6	aqu
tsk_7	aqu
tsk_8	aqu
tsk_9	aqu
tsk_10	pol
tsk_11	pol
tsk_12	pol
tsk_13	pol
tsk_14	pol
tsk_15	pol
tsk_16	pol
tsk_17	pol
tsk_18	pol
tsk_19	pol

###
### 2. Compute PIXY ------
###

### For Swiss aqu, Swiss pol (note that the empirical PIXY estimates are calculated or FINNISH aqu, Swiss pol; modelled demographic history not available for Finnish aqu / Swiss pol pair), 
### pixy_sim.sh

#!/bin/bash -l
#SBATCH -J pixy_sim
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/msprime/pixy/logs/pixy_sim.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/msprime/pixy/logs/pixy_sim.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

cd /scratch/project_2001443/barriers_introgr_formica/msprime/pixy/

pixy --stats pi dxy fst \
--vcf /scratch/project_2001443/barriers_introgr_formica/msprime/chr_test/sim_1/output.vcf.gz \
--populations groupfiles/group_pixy_sim.tab \
--window_size 10000 \
--output_prefix aqu_swi_pol_swi_sim \
--n_cores 8
