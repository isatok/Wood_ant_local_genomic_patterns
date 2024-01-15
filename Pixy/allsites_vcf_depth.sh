#This is to compute depth in 20kb windows from the allsites vcf to control for collapsed paralogues when detecting co-elevated Fst/Dxy regions

cd /scratch/project_2001443/barriers_introgr_formica/pixy

REFPATH=/scratch/project_2001443/reference_genome
VCF=/scratch/project_2001443/barriers_introgr_formica/gvcf/allsamples_filtered.vcf.gz
#generate regions with script from https://github.com/freebayes/freebayes/blob/master/scripts/fasta_generate_regions.py

## 2. Split ref in 20 kb regions
sinteractive --account project_2001443 --mem 2000
module load freebayes #v. 1.3.6
fasta_generate_regions.py $REFPATH/Formica_hybrid_v1_wFhyb_Sapis.fa.fai 20000 > Formica_hybrid_v1_20kb_regions.tmp

#remove regions that are not used: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00 
grep -v -e mtDNA -e wFhyb -e Spiroplasma Formica_hybrid_v1_20kb_regions.tmp > Formica_hybrid_v1_20kb_regions.txt ; rm Formica_hybrid_v1_20kb_regions.tmp

#subset VCF such that have no header, and only 1st (CHROM), 2nd (POS), and 8th (INFO) fields; from INFO only the number following DP 
bcftools query -f'[%CHROM:%POS %INFO/DP]' $VCF > allsamples_filtered_DP.txt

#get these info
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/allsamples_filtered_DP.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy/depth'

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/Formica_hybrid_v1_20kb_regions.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy/depth'





