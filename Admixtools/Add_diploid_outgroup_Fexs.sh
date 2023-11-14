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

# Combine outgroup + sample VCFs
VCFOUT=/scratch/project_2001443/barriers_introgr_formica/admixtools/all_samples_DP8_wFexs.vcf.gz

bcftools index -t $vcfexs

bcftools merge -Ou \
  $VCFIN \
  $vcfexs -Oz -o $VCFOUT

bcftools index -t $VCFOUT
