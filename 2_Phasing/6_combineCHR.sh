## Combine all chromosomes
cd $SCRATCH/vcf/phasing/shapeit
ls *vcf.gz > shapeit.file.list
bcftools concat -Oz -f shapeit.file.list -o all_samples.minDP8.AN10percMiss.mac2.whap.shapeit.allScafs.vcf.gz

# Remove temp files
rm -rf *it.Scaffold*vcf*
rm -rf ../whatshap/*vcf*

### 6_combineCHR.sh (END)
