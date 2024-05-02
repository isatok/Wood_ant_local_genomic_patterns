#Create needed directories:

mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/
mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/
mkdir /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/

cd /scratch/project_2001443/barriers_introgr_formica/gvcf/

#make a gvcf.bam.list
#find . -type f -name "*.bam" -exec realpath {} \; > gvcf.bam.list

#make a gvcf.bam.all.list for snp calling with all samples. filter out samples with poor data only after snp calling but before per-site 50% missing data filter
cd /scratch/project_2001443/barriers_introgr_formica/bam/bam_all/
find . -type f -name "*.bam" -exec realpath {} \; > /scratch/project_2001443/barriers_introgr_formica/gvcf/gvcf.bam.all.list 
#Here only RN417 and 121 are excluded as they were not in the SNP calling anyway due to bad sequencing data quality


### Batch script to make a gvcf: SNP calling and filtering ### 


#!/bin/bash -l
#SBATCH -J allsites_vcf_highcovsamples_mac2
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/allsites_vcf_highcovsamples_mac2_%j_%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/allsites_vcf_highcovsamples_mac2_%j_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 72:00:00
#SBATCH -p small
#SBATCH --array=1-26
#SBATCH --ntasks 1
#SBATCH --mem=8G
#SBATCH --mail-type=END

#
# 0. Load modules and set directory and variables -----------------------------
#

module load biokit

# go to directory where raw gvcf for each scaffold will be produced
cd /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/

# Get scaffold ID
REF=/scratch/project_2001443/reference_genome
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${REF}/scaffold_no0300.list)

# write gvcf
echo "Writing gvcf for ${scaffold} ..."

bcftools mpileup -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -b /scratch/project_2001443/barriers_introgr_formica/gvcf/gvcf.bam.all.list \
  -r ${scaffold} | bcftools call -m -Oz -f GQ -o /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/${scaffold}_allsamples.vcf.gz

# Set DP thresholds & filter missing data
#Max depth threshold: Calculate average of the mean depths and multiply by two -> follow the logic from the variant site vcf filtering where 
#sites with >2x mean ind depth are set as missing.
# IDEPTHPATH=/scratch/project_2001443/barriers_introgr_formica/vcf/filt
# cat $IDEPTHPATH/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth \
# | grep -v "INDV" | awk '{ total += $3; count++ } END { print total/count }' #15.6196
# 15.6196*2= ca. 31,24; set max avg depth to 31x.

mindp=206 #2*103
maxdp=3193 #31*103

echo "Filtering gvcf ..."
bcftools view -Ou ${scaffold}_allsamples.vcf.gz | bcftools filter -e "DP < ${mindp} | DP > ${maxdp}" --set-GTs . -Oz > ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz &&
bcftools index -t ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz

echo "Removing low coverage individuals ..."
vcftools --gzvcf ${scaffold}_allsamples_meanDP${mindp}-${maxdp}.vcf.gz --remove-indv 110-FaquH  --remove-indv RN418 --remove-indv 105-FaquH --remove-indv 54-Frufa \
--remove-indv s353 --remove-indv s354 --remove-indv RN421 --remove-indv RN425 --remove-indv RN426 --remove-indv RN422 \
--recode --recode-INFO-all --stdout | bgzip > ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz
bcftools index -t ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz

echo "Filtering missing data per site ..."
bcftools view -Ou ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}.vcf.gz | bcftools filter -e "F_MISSING > 0.5" -Oz > \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz


# two vcfs will be made, one with only invariant sites and one with only variants\
    # these two are then combine to make a vcf will all sites (gvcf)
# to combine they must both be indexed then combined with bcftools concat

echo "Creating vcf for invariant sites ..."
# Create a filtered VCF containing only invariant sites
vcftools --gzvcf ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz \
--max-maf 0 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz

echo "Creating vcf for variant sites ..."
# Create a filtered VCF containing only variant sites, keep SNPqual >= 30 and filter out HW excess and singletons
vcftools --gzvcf ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz --mac 2 --hwe 0.001 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

# Index both vcfs using tabix
tabix ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz
tabix ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

echo "Combining vcfs for ${scaffold} ..."
# Combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz \
${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz \
-Oz > ${scaffold}_highcovsamples_filtered.vcf.gz

# Index for concat afterwards
tabix ${scaffold}_highcovsamples_filtered.vcf.gz

echo "Removing temp files ..."
# Remove temp files
rm ${scaffold}_allsamples_meanDP${mindp}-${maxdp}_maxNA50perc*
rm ${scaffold}_allsamples_meanDP${mindp}-${maxdp}*
rm ${scaffold}_highcovsamples_meanDP${mindp}-${maxdp}*

# End bash script. ----------------------------





### Interactive script to combine the per scaffold gvcf files ### ----------------------------


# Combine VCFs
cd /scratch/project_2001443/barriers_introgr_formica/gvcf/raw/

#ls *allsamples_filtered.vcf.gz > all_filtered_vcfs.list
ls *_highcovsamples_filtered.vcf.gz > highcovsamples_filtered_vcfs.list

sinteractive...
module load biokit

bcftools concat \
--allow-overlaps \
-f highcovsamples_filtered_vcfs.list \
-Oz > ../highcovsamples_filtered.vcf.gz

tabix highcovsamples_filtered.vcf.gz


### Filter out duplicate sites from the all-sites vcf highcovsamples_filtered.vcf.gz ### ----------------------------


### remove_duplicates_highcov.sh ### 
#!/bin/bash -l
#SBATCH -J remove_duplicates_highcov
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates_highcov.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/gvcf/logs/remove_duplicates_highcov.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001443/barriers_introgr_formica/gvcf

VCFIN=/scratch/project_2001443/barriers_introgr_formica/gvcf/highcovsamples_filtered.vcf.gz

echo "counting the initial number of sites in highcovsamples_filtered.vcf.gz ..."
bcftools index -n $VCFIN

#Filter out duplicate sites (achieve this by filtering out all that are mnp's and not homozygous "ref" or "snp" type)

echo "counting the amount of duplicates in highcovsamples_filtered.vcf.gz (divide this number by two)..."
echo $VCFIN
gunzip -c $VCFIN | cut -f2 | uniq -D | wc -l

#Keep only sites hom for reference allele or snps
bcftools filter --threads 8 -Oz -e 'TYPE!="ref" && TYPE!="snp"' -m+ $VCFIN > highcovsamples_filtered_noindels.vcf.gz
bcftools index -t highcovsamples_filtered_noindels.vcf.gz

echo "counting how many sites remain after TYPEref and TYPEsnp filtering..."
bcftools index -n highcovsamples_filtered_noindels.vcf.gz

echo "counting how many duplicates are left (divide the number by two)..."
gunzip -c highcovsamples_filtered_noindels.vcf.gz | cut -f2 | uniq -D | wc -l

echo "checking examples of what kind of variants they are - mnps? look at these coordinates in the highcovsamples_filtered_noindels.vcf.gz file..."
gunzip -c highcovsamples_filtered_noindels.vcf.gz | grep "Scaffold01" | cut -f2 | uniq -D | head     # e.g. xx, yy, zz
#bcftools view -H -r Scaffold01:5047119 highcovsamples_filtered_noindels.vcf.gz #They are mnps

#Always the first one is the snp and the second one mnp. Remove all second instances
echo "removing all second instances of duplicate sites (expecting the first ones to be the real snps)..."
bcftools norm --rm-dup all --threads 8 -Oz highcovsamples_filtered_noindels.vcf.gz > highcovsamples_filtered_noindels_rmdup.vcf.gz
bcftools index -t highcovsamples_filtered_noindels_rmdup.vcf.gz

echo "counting how many sites remain after duplicate filtering (add to this the initial #of duplicates - does it match with the input vcf site count?) ..."
bcftools index -n highcovsamples_filtered_noindels_rmdup.vcf.gz 

echo "counting how many duplicates remain (divide by two)..."
gunzip -c highcovsamples_filtered_noindels_rmdup.vcf.gz | cut -f2 | uniq -D | wc -l

### END.
