###
### Prep input file w/ outgroup -------------------------------------------------------------------
###

# Create one BED file containing the positions of all phased SNPs
# This BED file will be used as an argument to call genotypes only at these positions

cd /scratch/project_2001443/vcf/phasing/shapeit

PHASEDVCF=/scratch/project_2001443/vcf/phasing/shapeit/all_samples.minDP8.AN10percMiss.mac2.whap.shapeit.allScafs.vcf.gz

# Extract first 2 columns (chr - pos)
gunzip -c $PHASEDVCF | grep -v '^#' | cut -f1,2 > phased_sites.tmp

# Duplicate pos column for BED format
cut -f2 phased_sites.tmp > tmp
paste -d '\t' phased_sites.tmp tmp > tmp2

# Compute BED range (BED is 0-based)
awk -F"  " 'BEGIN{OFS="  "} {$2=$2-1;print $0}' tmp2 > phased_sites.bed
rm phased_sites.tmp tmp tmp2


###
### SNP calling for F. exsecta --------------------------------------------------------------------
###

# Instead of Freebayes, it is simpler / quicker to genotype using bcftools
# Run on an interactive node; check that ploidy 1 flag works as expected (if a site is polymorphic, should get the major allele)

sinteractive
cd /scratch/project_2001443/vcf/phasing/shapeit

module load biokit

bam=/scratch/project_2001443/Fexsecta_bam/Fexs_nodupl.bam
ref=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa
phbed=/scratch/project_2001443/vcf/phasing/shapeit/phased_sites.bed
#bcftools mpileup -f $ref -R $phbed $bam | bcftools call --ploidy 1 -m -Oz -f GQ -o Fexs_nodupl_phased_sites_tmp.vcf.gz #now this may be overcommented by mistake? should indeed get a haplodised outgroup
bcftools mpileup -f $ref -R $phbed $bam | bcftools call -m -Oz -f GQ -o Fexs_nodupl_phased_sites_tmp.vcf.gz

# Edit header (instead of adding the read group, here the name of the sample will be "Fexs")
echo "Fexs" >> out.name
bcftools reheader -s out.name Fexs_nodupl_phased_sites_tmp.vcf.gz -o Fexs_nodupl_phased_sites.vcf.gz
rm out.name Fexs_nodupl_phased_sites_tmp.vcf.gz

# Combine outgroup + sample VCFs
PHASEDVCF=/scratch/project_2001443/vcf/phasing/shapeit/all_samples.minDP8.AN10percMiss.mac2.whap.shapeit.allScafs.vcf.gz

bcftools merge -Ou \
  $PHASEDVCF \
  Fexs_nodupl_phased_sites.vcf.gz -Oz -o phased_with_outgroup.vcf.gz

  #Set missing genotypes (which are only in exsecta) to . instead of ./.
zcat phased_with_outgroup.vcf.gz | awk '{gsub(/\.\/\./,"."); print $0}' | bgzip > phased_with_outgroup_gtfix.vcf.gz
