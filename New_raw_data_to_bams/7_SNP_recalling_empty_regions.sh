#Inspect the called VCF for SNP-empty 50kb regions (due to time-out) and do more calling in these regions after normalization


sinteractive --account project_2001443 --mem 6000

cd /scratch/project_2001443/barriers_introgr_formica/vcf

# the original SNP-calling regions file: /scratch/project_2001443/vcf/allas_vcf_check/filt/Formica_hybrid_v1_50kb_regions.bed
# the SNP-empty-regions-file for Satokangas et al. (2023): /scratch/project_2001443/vcf/allas_vcf_check/filt/regions_wout_SNPs_INA.bed 

module load biokit

#index the vcf's & check how many sites (check file name)
bcftools index -t all_samples_XXX.normalized.vcf.gz
bcftools index -n all_samples_XXX.normalized.vcf.gz

####CONTINUE EDITING FROM HERE ######

# DONE # 3. make a bed file out of the .txt region file simply by renaming it, and then replace : & / with tabs.
#(in sed, "g" stands for global replacement i.e. all occurrences):
cp /scratch/project_2001443/reference_genome/Formica_hybrid_v1_50kb_regions.txt /scratch/project_2001443/vcf/filt/Formica_hybrid_v1_50kb_regions.bed
cp /scratch/project_2001443/reference_genome/regions_wout_SNPs_INA.txt /scratch/project_2001443/vcf/filt/regions_wout_SNPs_INA.bed

sed 's/:/\t/g' Formica_hybrid_v1_50kb_regions.bed | sed 's/-/\t/g' > Formica_hybrid_v1_50kb_regions_t.bed
sed 's/:/\t/g' regions_wout_SNPs_INA.bed | sed 's/-/\t/g' > regions_wout_SNPs_INA_t.bed

# 4. Count per window for newly called "SNP-empty-regions" (2) file:
#trying to do this as a batch job, as sinteractive would need more than the automatic resources, and don't know how much
bedtools coverage -a regions_wout_SNPs_INA_t.bed -b all_samples_2.normalized.vcf.gz -counts > 50kb_regions_SNP_counts_2.tab

# 5. merge the old SNP-call-file and the newer empty-regions-SNP-call-file
bcftools merge --force-samples all_samples.normalized.vcf.gz all_samples_2.normalized.vcf.gz -Oz -o all_samples_3.normalized.vcf.gz

# 6. index the new merged file (see step 2.)

# 7. Count SNPs per 50kb window for merged (3) file:
bedtools coverage -a 50kb_regions.bed -b all_samples_3.normalized.vcf.gz -counts > 50kb_regions_SNP_counts_3.tab

# 8. Plot the counts along the genome, check especially region number 130, as it was cancelled due to timeout

# 9. We'll continue with filtering later on, when Pierre has done the pipeline.
