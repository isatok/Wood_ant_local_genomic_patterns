## Combine all chromosomes
module load biokit

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
ls *vcf.gz > shapeit.file.list
bcftools concat -Oz -f shapeit.file.list -o all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.shapeit.allScafs.vcf.gz

# Remove temp files
rm -rf *it.Scaffold*vcf*
rm -rf ../whatshap/*vcf*

### 6_combineCHR.sh (END)
