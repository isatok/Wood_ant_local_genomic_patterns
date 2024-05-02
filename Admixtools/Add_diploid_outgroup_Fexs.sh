
### Not used for the current manuscript (updated in May 2024) ###


###
### Prep input file w/ outgroup -------------------------------------------------------------------
###

# Create one BED file containing the positions of all SNPs
# This BED file will be used as an argument to call genotypes only at these positions

cd /scratch/project_2001443/barriers_introgr_formica/admixtools

#A VCF containing all old and new sequenced samples for which we have even some ok data (ie excl. RN117(?))
VCFIN=/scratch/project_2001443/barriers_introgr_formica/vcf/filt/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.vcf.gz

# Extract first 2 columns (chr - pos)
gunzip -c $VCFIN | grep -v '^#' | cut -f1,2 > vcf_sites.tmp

# Duplicate pos column for BED format
cut -f2 vcf_sites.tmp > tmp
paste -d '\t' vcf_sites.tmp tmp > tmp2

# Compute BED range (BED is 0-based); to insert tabs (in between ""; 2x), press [ctrl] + v + [tab]
awk -F"  " 'BEGIN{OFS="  "} {$2=$2-1;print $0}' tmp2 > vcf_sites.bed
rm vcf_sites.tmp tmp tmp2


###
### SNP calling for F. exsecta --------------------------------------------------------------------
###

# Instead of Freebayes, it is simpler / quicker to genotype using bcftools
# Run on an interactive node; (note: bcftools call --ploidy 1 retains info about both ref and alt alleles, but indicates in the GT field that it's haploid & correctly chooses the major allele)

bam=/scratch/project_2001443/Fexsecta_bam/Fexs_nodupl.bam
ref=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa
phbed=/scratch/project_2001443/barriers_introgr_formica/admixtools/vcf_sites.bed
vcfexs=Fexs_nodupl_sites.vcf.gz

bcftools mpileup -f $ref -R $phbed $bam | bcftools call --ploidy 1 -m -Oz -f GQ -o Fexs_nodupl_sites_tmp.vcf.gz #took ca. 15 minutes

# Edit header (instead of adding the read group, here the name of the sample will be "Fexs")
echo "Fexs" >> out.name
bcftools reheader -s out.name Fexs_nodupl_sites_tmp.vcf.gz -o $vcfexs

rm out.name Fexs_nodupl_sites_tmp.vcf.gz

###REMOVE EXSECTA DUPLICATE SITES, START FROM HERE###

#How many sites to begin with?
bcftools index -n $VCFIN #708783 in the original vcf; also wc -l $phbed 708783 is in line
bcftools index -n $vcfexs #696424 in the exsecta vcf

# Filter out duplicate sites (achieve this by filtering out all that are mnp's and not homozygous "ref" or "snp" type)

gunzip -c $vcfexs | cut -f2 | uniq -D | wc -l #Altogether 1814/2 = 907 duplicates in the unfiltered exsecta file

#Keep only sites hom for reference allele or snps
bcftools filter --threads 8 -Oz -e 'TYPE!="ref" && TYPE!="snp"' -m+ $vcfexs > Fexs_nodupl_sites_noindels.vcf.gz
bcftools index -t Fexs_nodupl_sites_noindels.vcf.gz 
bcftools index -n Fexs_nodupl_sites_noindels.vcf.gz  #695657 -> 767 sites removed

gunzip -c Fexs_nodupl_sites_noindels.vcf.gz | cut -f2 | uniq -D | wc -l #280 /2 =140 #We still have 140 duplicates
gunzip -c Fexs_nodupl_sites_noindels.vcf.gz | grep "Scaffold01" | cut -f2 | uniq -D     # e.g. 2442227, 5047119, 5422749
bcftools view -H -r Scaffold01:5047119 Fexs_nodupl_sites_noindels.vcf.gz #They are mnps

#Always the first one is the snp and the second one mnp. Remove all second instances
bcftools norm --rm-dup all --threads 8 -Oz Fexs_nodupl_sites_noindels.vcf.gz > Fexs_nodupl_sites_noindels_rmdup.vcf.gz
bcftools index -t Fexs_nodupl_sites_noindels_rmdup.vcf.gz
bcftools index -n Fexs_nodupl_sites_noindels_rmdup.vcf.gz #695517 snps remain
gunzip -c Fexs_nodupl_sites_noindels_rmdup.vcf.gz | cut -f2 | uniq -D | wc -l #no duplicates left; 695517(remaining)+907(duplicates)=696424 sites, equals to original amount of sites in $vcfex

#Merge the actual dataset and exsecta vcfs

VCFOUT=/scratch/project_2001443/barriers_introgr_formica/admixtools/all_samples_DP8_wFexs_dedupl.vcf.gz

#bcftools index -t $vcfexs

bcftools merge -Ou \
  $VCFIN \
  Fexs_nodupl_sites_noindels_rmdup.vcf.gz -Oz -o $VCFOUT

bcftools index -t $VCFOUT
bcftools index -n $VCFOUT # 708783 variants; all good
bcftools query -l $VCFOUT | wc -l # 104 samples -> remember to remove low-q inds

# Change GTs in haploid Fexs sample from ./. into . ### NOTE THAT IF THE NUMBER OF SAMPLES CHANGES THIS NEEDS TO CHANGE. NOW ALTOGETHER 113 COLS (9 + 104); FEXS IS THE LAST ONE. ###
zcat all_samples_DP8_wFexs_dedupl.vcf.gz | awk 'BEGIN {OFS=FS="\t"} {gsub(/\.\/\./,".", $113); print $0}' | bgzip >  all_samples_DP8_wFexs_dedupl_gtfix.vcf.gz
bcftools index -t all_samples_DP8_wFexs_dedupl_gtfix.vcf.gz

# The vcf file originating from all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.vcf.gz., i.e. unphased and
# no mac-filtering, but outgroup merged and duplicates taken care of, is
# all_samples_DP8_wFexs_dedupl_gtfix.vcf.gz at /scratch/project_2001443/barriers_introgr_formica/admixtools. This should be used for ADMIXTOOLS and Fbranch.
