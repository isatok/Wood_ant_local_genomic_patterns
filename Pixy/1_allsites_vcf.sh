#Create needed directories:
mkdir /scratch/project_2001443/gvfc/ ### NOT LIKE THIS. REMOVE IF EXISTS
mkdir /scratch/project_2001443/gvfc/logs/ ### NOT LIKE THIS. REMOVE IF EXISTS

mkdir /scratch/project_2001443/barriers_introgression_formica/gvfc/
mkdir /scratch/project_2001443/barriers_introgression_formica/gvfc/raw/
mkdir /scratch/project_2001443/barriers_introgression_formica/gvfc/logs/


#Make or copy needed lists:
scaffold.list
sample.list (no #110, (#54,) #105, #RN417) (?) Currently IDEPTH calculated for 101 samples !!
bam.list



#!/bin/bash -l
#SBATCH -J allsites_vcf
#SBATCH -o /scratch/project_2001443/gvfc/logs/allsites_vcf_%j_%a.out
#SBATCH -e /scratch/project_2001443/gvfc/logs/allsites_vcf_%j_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 72:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 1
#SBATCH --mem=8G


#
# 0. Load modules and set directory and variables -----------------------------
#

module load biokit

# go to directory with bam files in 
cd /scratch/project_2001443/barriers_introgr_formica/xxx ######CHECK

# Get scaffold ID
REF=/scratch/project_2001443/reference_genome ######CHECK
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p $REF/scaffold.list) ######CHECK

# write gvcf
echo Writing gvcf for $scaffold ...

bcftools mpileup -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -b /scratch/project_2001443/bam/mkduple_RG_clip_bams_246samples_onlymales.list  \ ##BAMLIST HERE
  -r ${scaffold} | bcftools call -m -Oz -f GQ -o /scratch/project_2001443/barriers_introgr_formica/gvfc/raw/${scaffold}_allsamples.vcf.gz


# go to directory with gvcf for each scaffold
cd /scratch/project_2001443/barriers_introgression_formica/gvfc/raw/

# Set DP thresholds & filter missing data
#Max depth threshold: Calculate average of the mean depths and multiply by two -> follow the logic from the variant site vcf filtering where sites with >2x mean ind depth are set as missing. Bc of big variance add some?
# IDEPTHPATH=/scratch/project_2001443/barriers_introgr_formica/vcf/filt
# awk '{ total += $3; count++ } END { print total/count }' $IDEPTHPATH/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth #15.4694
# 15.4694*2= ca. 31; set max avg depth to 40x?
mindp= # 2*101
maxdp= # 40*101

echo Filtering gvcf ....
bcftools view -Ou ${scaffold}_all350samples_noDiploidMales.vcf.gz | bcftools filter -e "DP < $mindp  | DP < $maxdp" --set-GTs . -Ou | bcftools filter -e "F_MISSING > 0.5" > ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz

# two vcfs will be made, one with only invariant sites and one with only variants\
    # these two are then combine to make a vcf will all sites (gvcf)
# to combine they must both be indexed then combined with bcftools concat

echo Creating vcf for invariant sites ...
# Create a filtered VCF containing only invariant sites
vcftools --gzvcf ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz \
--max-maf 0 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz

echo Creating vcf for variant sites ...
# Create a filtered VCF containing only variant sites, keep SNPqual >= 30 and filter out HW excess
vcftools --gzvcf ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz --mac 1 --hwe 0.001 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz



# Index both vcfs using tabix
tabix ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz
tabix ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

echo Combining vcfs for $scaffold ...
# Combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz \
${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz \
-Oz > ${scaffold}_all350samples_noDiploidMales_filtered.vcf.gz

# Index for concat afterwards
tabix ${scaffold}_all350samples_noDiploidMales_filtered.vcf.gz

echo Removing temp files ...
# Remove temp files
# rm ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc*
rm ${scaffold}_all350samples_noDiploidMales_meanDP${mindp}-${maxdp}_maxNA50perc*

# End bash script
