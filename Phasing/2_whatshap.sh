#!/bin/bash -l
#SBATCH -J whap
#SBATCH -o /scratch/project_2001443/vcf/phasing/whatshap/logs/whatshap_%a.out
#SBATCH -e /scratch/project_2001443/vcf/phasing/whatshap/logs/whatshap_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-89
#SBATCH --ntasks 1
#SBATCH --mem=3G

# conda mambo
module load biokit
module load bioconda/3
conda activate my_seqdata

cd $SCRATCH/vcf/phasing

# Define paths
BAMDIR=$SCRATCH/bam/nodupl_RG_clip/renamed
masterVCF=$SCRATCH/vcf/filt/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.mac2.vcf.gz
REF=$SCRATCH/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa

# Define focal ind
ind=$(sed -n "$SLURM_ARRAY_TASK_ID"p ind.list)

# 1. Extract single individual from VCF & fix floats in GQ field (or issue later when parsing with whatshap)
bcftools view -s $ind $masterVCF -Ov | perl -npe 's/(.\/.)\:.+\:(.+\:.+\:.+\:.+\:.+\:.+\:.+)/$1\:99\:$2/' | grep -v '^Scaffold00' | bgzip > whatshap/${ind}.GQfixed.vcf.gz &&\
bcftools index -t whatshap/${ind}.GQfixed.vcf.gz

# 2. Phase
whatshap phase -o whatshap/${ind}.phased.vcf.gz --reference $REF \
               whatshap/${ind}.GQfixed.vcf.gz $BAMDIR/${ind}_nodupl_*.bam &&\
whatshap stats whatshap/${ind}.phased.vcf.gz --tsv=whatshap/${ind}.phased.tsv

### end
2_whatshap.sh (END)
