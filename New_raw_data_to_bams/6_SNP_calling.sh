# 6_SNP_calling.sh


cd /scratch/project_2001443/vcf


###
### Copy to the major bam folder a) F. paralugubris bam files and b) bam files from the previous project
###

a) F. paralugubris bam files
/scratch/project_2001443/paralugubris/bam/nodupl_RG_clip
b) bam files from the previous project
/scratch/project_2001443/barriers_introgr_formica/bam_semipermeable

###
### Create additional input files
###

# 1. BAM list w/out samples 121 and RNXXX (excluded due to low quality)
ls /scratch/project_2001443/bam/nodupl_RG_clip/*.bam > bam.tmp
grep -v -e 121 -e RNXXX bam.tmp > bam.list ; rm bam.tmp

# 2. Split ref in 50 kb regions to speed up the analysis
module load freebayes
fasta_generate_regions.py Formica_hybrid_v1.fa.fai 50000 > Formica_hybrid_v1_50kb_regions.tmp




###
### Full SNP calling
###

# --skip-coverage = total DP combined across all samples; assuming max 400x per sample per region, 400X per sample * 101 samples = 40400X
# use screen for the SNP calling https://linuxize.com/post/how-to-use-linux-screen/?utm_content=cmp-true

module load freebayes # version 2023: v1.3.6

cd /scratch/project_2001443/vcf
REF=/scratch/project_2001443/reference_genome
RES=/scratch/project_2001443/vcf

/appl/soft/bio/bioconda/miniconda3/envs/freebayes/bin/freebayes-puhti \
  -time 72 \
  -regions $REF/Formica_hybrid_v1_50kb_regions.txt \
  -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L $RES/bam.list \
  -k --genotype-qualities --skip-coverage 40400 \
  -out $RES/raw/all_samples_raw.vcf




###
### Compress & sort
###

#!/bin/bash -l
#SBATCH -J comp_sort
#SBATCH -o /scratch/project_2001443/vcf/raw/compress_sort.out
#SBATCH -e /scratch/project_2001443/vcf/raw/compress_sort.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001443/vcf/raw
# Compress
/appl/soft/bio/samtools/htslib-1.9/bgzip all_samples_raw.vcf

# Sort
bcftools sort -m 1G -O z -o all_samples.vcf.gz -T ./tmp_sort all_samples_raw.vcf.gz

# Index
/appl/soft/bio/samtools/htslib-1.9/tabix -p vcf all_samples.vcf.gz
bcftools index -n all_samples.vcf.gz

gunzip -c all_samples.vcf.gz | grep -v "^#" | cut -f 1 | uniq -c > site_count_per_chromosome.tab




### END
