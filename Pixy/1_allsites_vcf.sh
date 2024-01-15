#Create needed directories:

mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/
mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/
mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/

#make a gvcf.bam.list
find . -type f -name "*.bam" -exec realpath {} \; > gvcf.bam.list


### Batch script to make a gvcf: SNP calling and filtering ### 


#!/bin/bash -l
#SBATCH -J allsites_vcf
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/allsites_vcf_%j_%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/allsites_vcf_%j_%a.err
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

# go to directory where raw gvcf for each scaffold will be produced
cd /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/

# Get scaffold ID
REF=/scratch/project_2001443/reference_genome
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${REF}/scaffold.list)

# write gvcf
echo "Writing gvcf for ${scaffold} ..."

bcftools mpileup -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -b /scratch/project_2001443/barriers_introgr_formica/gvcf/gvcf.bam.list
  -r ${scaffold} | bcftools call -m -Oz -f GQ -o /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/${scaffold}_allsamples.vcf.gz

# Set DP thresholds & filter missing data
#Max depth threshold: Calculate average of the mean depths and multiply by two -> follow the logic from the variant site vcf filtering where 
#sites with >2x mean ind depth are set as missing. Bc of big variance add some?
# IDEPTHPATH=/scratch/project_2001443/barriers_introgr_formica/vcf/filt
# cat $IDEPTHPATH/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth | grep -v "105-FaquH" | grep -v "110-FaquH" \
# | grep -v "INDV" | awk '{ total += $3; count++ } END { print total/count }' #15.5604
# 15.5604*2= ca. 31,12; set max avg depth to 32x.

mindp=202 #2*101
maxdp=3232 #32*101

echo "Filtering gvcf ..."
bcftools view -Ou ${scaffold}_allsamples.vcf.gz | bcftools filter -e "DP < ${mindp} | DP > ${maxdp}" --set-GTs . -Ou | bcftools filter -e "F_MISSING > 0.5" -Oz > \
${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz

# two vcfs will be made, one with only invariant sites and one with only variants\
    # these two are then combine to make a vcf will all sites (gvcf)
# to combine they must both be indexed then combined with bcftools concat

echo "Creating vcf for invariant sites ..."
# Create a filtered VCF containing only invariant sites
vcftools --gzvcf ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz \
--max-maf 0 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz

echo "Creating vcf for variant sites ..."
# Create a filtered VCF containing only variant sites, keep SNPqual >= 30 and filter out HW excess
vcftools --gzvcf ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz --mac 1 --hwe 0.001 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > \
${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

# Index both vcfs using tabix
tabix ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz
tabix ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

echo "Combining vcfs for ${scaffold} ..."
# Combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz \
${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz \
-Oz > ${scaffold}_allsamples_filtered.vcf.gz

# Index for concat afterwards
tabix ${scaffold}_allsamples_filtered.vcf.gz

echo "Removing temp files ..."
# Remove temp files
rm ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc*

# End bash script. ----------------------------




### Interactive script to combine the per scaffold gvcf files ### 


# Combine VCFs
cd /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/
#for file in *allsamples_filtered.vcf.gz ; do tabix $file ; done

ls *allsamples_filtered.vcf.gz > all_filtered_vcfs.list

sinteractive...
module load biokit

bcftools concat \
--allow-overlaps \
-f all_filtered_vcfs.list \
-Oz > ../allsamples_filtered.vcf.gz

tabix allsamples_filtered.vcf.gz



### Filter out duplicate sites from the all-sites vcf ### TO BE DONE
### FOLLOW THIS PIPELINE MADE FOR EXSECTA (MODIFY AS NEEDED) ###

VCFIN=/scratch/project_2001443/barriers_introgr_formica/gvcf/allsamples_filtered.vcf.gz

#How many sites to begin with?
bcftools index -n $VCFIN #208.029.582

# Filter out duplicate sites (achieve this by filtering out all that are mnp's and not homozygous "ref" or "snp" type)

gunzip -c $VCFIN | cut -f2 | uniq -D | wc -l #Altogether 268.064 /2 = 134.032 duplicates in the unfiltered file

#Keep only sites hom for reference allele or snps
bcftools filter --threads 8 -Oz -e 'TYPE!="ref" && TYPE!="snp"' -m+ $VCFIN > allsamples_filtered_noindels.vcf.gz
bcftools index -t allsamples_filtered_noindels.vcf.gz
bcftools index -n allsamples_filtered_noindels.vcf.gz  #xx -> yy sites removed

#gunzip -c allsamples_filtered_noindels.vcf.gz | cut -f2 | uniq -D | wc -l #xx /2 =yy #We still have yy duplicates
#gunzip -c allsamples_filtered_noindels.vcf.gz | grep "Scaffold01" | cut -f2 | uniq -D     # e.g. xx, yy, zz
#bcftools view -H -r Scaffold01:5047119 allsamples_filtered_noindels.vcf.gz #They are mnps

#Always the first one is the snp and the second one mnp. Remove all second instances
bcftools norm --rm-dup all --threads 8 -Oz allsamples_filtered_noindels.vcf.gz > allsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -t allsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -n allsamples_filtered_noindels_rmdup.vcf.gz #xx snps remain
gunzip -c allsamples_filtered_noindels_rmdup.vcf.gz | cut -f2 | uniq -D | wc -l #no duplicates left; xx(remaining)+yy(duplicates)=zz sites, equals to original amount of sites in $vcfex



### Try to do this as a batch job as the invariant file is so huge, if needed; remove_duplicates.sh ### MODIFY THIS BATCH JOB IF NEEDED

#!/bin/bash -l
#SBATCH -J remove_duplicates
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates.err
#SBATCH --account=project_2001443
#SBATCH -t 72:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001443/barriers_introgr_formica/gvcf

VCFIN=/scratch/project_2001443/barriers_introgr_formica/gvcf/allsamples_filtered.vcf.gz

#How many sites to begin with?
bcftools index -n $VCFIN #208.029.582

# Filter out duplicate sites (achieve this by filtering out all that are mnp's and not homozygous "ref" or "snp" type)

gunzip -c $VCFIN | cut -f2 | uniq -D | wc -l #Altogether xx/2 = yy duplicates in the unfiltered exsecta file

#Keep only sites hom for reference allele or snps
bcftools filter --threads 8 -Oz -e 'TYPE!="ref" && TYPE!="snp"' -m+ $VCFIN > allsamples_filtered_noindels.vcf.gz
bcftools index -t allsamples_filtered_noindels.vcf.gz
bcftools index -n allsamples_filtered_noindels.vcf.gz  #xx -> yy sites removed

gunzip -c allsamples_filtered_noindels.vcf.gz | cut -f2 | uniq -D | wc -l #xx /2 =yy #We still have yy duplicates
gunzip -c allsamples_filtered_noindels.vcf.gz | grep "Scaffold01" | cut -f2 | uniq -D     # e.g. xx, yy, zz
bcftools view -H -r Scaffold01:5047119 allsamples_filtered_noindels.vcf.gz #They are mnps

#Always the first one is the snp and the second one mnp. Remove all second instances
bcftools norm --rm-dup all --threads 8 -Oz allsamples_filtered_noindels.vcf.gz > allsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -t allsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -n allsamples_filtered_noindels_rmdup.vcf.gz #xx snps remain
gunzip -c allsamples_filtered_noindels_rmdup.vcf.gz | cut -f2 | uniq -D | wc -l #no duplicates left; xx(remaining)+yy(duplicates)=zz sites, equals to original amount of sites in $vcfex

