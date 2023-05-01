sinteractive --account project_2001443 --mem 2000

#set vcf to be filtered for minor allele count
VCF=/scratch/project_2001443/vcf/filt/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.vcf.gz

#set minimum minor allele count to 2, i.e. remove singletons
vcftools --gzvcf $VCF --mac 2 --recode --stdout | gzip -c > out.mac2.vcf.gz

#bgzip the file and rename it properly
gunzip -c out.mac2.vcf.gz | bgzip > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.mac2.vcf.gz

#index the file
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.mac2.vcf.gz

#count the remaining non-singleton variants #1829565
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.mac2.vcf.gz
1_mac_filtering.sh (END)
