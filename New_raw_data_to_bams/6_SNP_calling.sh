# 6_SNP_calling.sh

cd /scratch/project_2001443/vcf


###
### Copy to the major bam folder a) F. paralugubris bam files and b) bam files from the previous project
###

FINALBAMDIR=/scratch/project_2001443/barriers_introgr_formica/bam_all

# a) F. paralugubris bam files (for Seifert)
cp /scratch/project_2001443/paralugubris/bam/nodupl_RG_clip/*.bam $FINALBAMDIR
cp /scratch/project_2001443/paralugubris/bam/nodupl_RG_clip/*.bam.bai $FINALBAMDIR

# b) bam files from the previous project (Satokangas et al 2023)
cp /scratch/project_2001443/barriers_introgr_formica/bam_semipermeable/*.bam $FINALBAMDIR
cp /scratch/project_2001443/barriers_introgr_formica/bam_semipermeable/*.bam.bai $FINALBAMDIR

# c) newly created bam files for this project
cp /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip/*.bam $FINALBAMDIR
cp /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip/*.bam.bai $FINALBAMDIR


###
### Create additional input files
###

# 1. BAM list w/out samples 121 and RN417 (RN417 excluded due to low quality; coverage only 2.16 after overlap correction)
ls $FINALBAMDIR/*.bam > $FINALBAMDIR/bam.tmp
grep -v -e 121 -e RN417 $FINALBAMDIR/bam.tmp > $FINALBAMDIR/bam.list ; rm $FINALBAMDIR/bam.tmp

## 2. Split ref in 50 kb regions to speed up the analysis #https://docs.csc.fi/apps/freebayes/
#sinteractive --account project_2001443 --mem 2000
#module load freebayes #v. 1.3.6 - different v. from earlier Satokangas et al 2023 pipeline; ok since now all data is re-prepared.
#cd /scratch/project_2001443/reference_genome
#fasta_generate_regions.py Formica_hybrid_v1_wFhyb_Sapis.fa.fai 50000 > Formica_hybrid_v1_50kb_regions.tmp

# 2. Split ref in equal data size regions to speed up the analysis #https://docs.csc.fi/apps/freebayes/

sinteractive --account project_2001443 --mem 2000
cd /scratch/project_2001443/barriers_introgr_formica/vcf
nano split_ref_by_bai_datasize.py #paste code from https://github.com/freebayes/freebayes/blob/master/scripts/split_ref_by_bai_datasize.py and change 'python3' into 'python'
module load biopythontools

python split_ref_by_bai_datasize.py \
-L /scratch/project_2001443/barriers_introgr_formica/bam_all/bam.list \
-r /scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa.fai \
-s 5000000 > /scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis_5e6_data_regions.fa.fai

#modify regions file so that Freebayes accepts it
awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' Formica_hybrid_v1_wFhyb_Sapis_5e6_data_regions.fa.fai > Formica_hybrid_v1_5e6_data_regions.tmp

#remove regions that require different SNP calling parameters if they are needed: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*)
cd /scratch/project_2001443/reference_genome/
#grep -v -e mtDNA -e wFhyb -e Spiroplasma Formica_hybrid_v1_50kb_regions.tmp > Formica_hybrid_v1_50kb_regions.txt ; rm Formica_hybrid_v1_50kb_regions.tmp
grep -v -e mtDNA -e wFhyb -e Spiroplasma Formica_hybrid_v1_5e6_data_regions.tmp > Formica_hybrid_v1_5e6_data_regions.txt ; rm Formica_hybrid_v1_5e6_data_regions.tmp


###
### Full SNP calling
###

# --skip-coverage = total DP combined across all samples; assuming max 400x per sample per region, 400X per sample * 103 samples = 41200X
# use screen for the SNP calling https://linuxize.com/post/how-to-use-linux-screen/?utm_content=cmp-true

screen -S snp_calling # create a named screen session
# Ctrl + a ? - list of commands
# Ctrl + a d - detach
screen -r # resume screen session


module load freebayes # version 2023: v1.3.6

cd /scratch/project_2001443/barriers_introgr_formica/vcf
REF=/scratch/project_2001443/reference_genome
RES=/scratch/project_2001443/barriers_introgr_formica/vcf
BAM=/scratch/project_2001443/barriers_introgr_formica/bam_all

freebayes-puhti \
  -time 72 \
  -mem 64 \
  -regions $REF/Formica_hybrid_v1_5e6_data_regions.txt \
  -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L $BAM/bam.list \
  -k --genotype-qualities --skip-coverage 41200 \
  -out $RES/raw/all_samples_raw.vcf

#MORE POTENTIAL OPTIONS IF NEEDED:

freebayes-puhti -regions regions.txt --limit-coverage 50 -k --genotype-qualities -f /scratch/project_2003480/patrick/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa  -L test_bams_last10_bam_bam22.list -out test_bams_last10_bam_bam22.vcf

--limit-coverage 50 #patrick has 5-10x, so 50x is repeat territory
-n best alleles
-E N=-1 disable complex variants - maybe use??


#########INA CONTINUE FROM HERE 3.8.2023###########


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
