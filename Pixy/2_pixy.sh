#### This script is to compute local Pi, Dxy and Fst estimates with PIXY ####

mkdir /scratch/project_2001443/barriers_introgr_formica/pixy
mkdir /scratch/project_2001443/barriers_introgr_formica/pixy/groupfiles
mkdir /scratch/project_2001443/barriers_introgr_formica/pixy/logs

cd /scratch/project_2001443/barriers_introgr_formica/pixy


###
### 1. Prepare group files for PIXY ------
###

# see at the bottom (takes a lot of space so instructions to make these are at the bottom of the script) #


###
### 2. Compute PIXY ------
###


###ONCE MORE 21.1.2024, with updated groupsfile to be similar w "simple" Fbranch & updated highcov vcf file###
###REMEMBER THAT W RUFA/POL COMPARISONS WE HAVE INTROGRESSION DATA FOR ONLY A SUBSET OF THE SAMPLES GROUPED HERE###

#/scratch/project_2001443/barriers_introgr_formica/pixy/groupfiles/group_pixy_100kb_5inds_highcovsamples.tab
#/scratch/project_2001443/barriers_introgr_formica/gvcf/highcovsamples_filtered_noindels_rmdup.vcf.gz

#!/bin/bash -l
#SBATCH -J pixy_100kb_highcov_5inds
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/pixy/logs/pixy_100kb_highcov_5inds.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/pixy/logs/pixy_100kb_highcov_5inds.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

cd /scratch/project_2001443/barriers_introgr_formica/pixy/

pixy --stats pi dxy fst \
--vcf /scratch/project_2001443/barriers_introgr_formica/gvcf/highcovsamples_filtered_noindels_rmdup.vcf.gz \
--populations groupfiles/group_pixy_100kb_5inds_highcovsamples.tab \
--window_size 100000 \
--output_prefix pixy_100kb_highcov_5inds \
--n_cores 8

#END.


#get the results for plotting
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/pixy_100kb_5inds*.txt' '/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy'

