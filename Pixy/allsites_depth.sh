### Compute bam coverage for the highcov individuals used in Jan 2024 ###

cd /scratch/project_2001443/bam/cov_highcovInds_5perSp

#make a list with filenames and bampaths

file=/scratch/project_2001443/bam/cov_highcovInds_5perSp/bam.highcovinds.list 

#This is the order of the individuals as how they are in the cleaned coverage file (according to the order in ind_cov.list created later)
#ind_cov_species.list:
chr chr
start start
end end
102-Frufa rufa
104-Faqu aqu
106-Faqu aqu
107-Faqu aqu
108-Flug lug
109-Faqu aqu
113-Flug lug
115-Flug lug
117-Fprat prat
11-Flug lug
120-Fprat prat
122-Fprat prat
123-Fprat prat
124-Flug lug
16-Frufa rufa
26-Frufa rufa
37-Frufa rufa
44-FaquH aqu
6-Fprat prat
72-Frufa rufa
CAGa_1w pol
CBCH1_1w pol
CBCH3_1w pol
NAZa_1w pol
VDa_1w pol



###
### Compute sequencing depth per sample (after duplicate removal) ------------------------------------------------------------------------
###



#!/bin/bash -l
#SBATCH -J mosdepth_highcovinds
#SBATCH -o /scratch/project_2001443/bam/cov_highcovInds_5perSp/logs/mosdepth_highcovinds%j.out
#SBATCH -e /scratch/project_2001443/bam/cov_highcovInds_5perSp/logs/mosdepth_highcovinds%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-25
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

BAMDIR=/scratch/project_2001443/barriers_introgr_formica/bam/bam_all
FINALDIR=/scratch/project_2001443/bam/cov_highcovInds_5perSp
FILELIST=/scratch/project_2001443/bam/cov_highcovInds_5perSp/bam.highcovinds.list
NAMELIST=/scratch/project_2001443/bam/cov_highcovInds_5perSp/nodupl_highcovinds_numonly.list

cd /scratch/project_2001443/bam/cov_highcovInds_5perSp

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILELIST)

# Get sample ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $NAMELIST)

mosdepth -t 1 -b 100000 -n $FINALDIR/${sample}_overlap_correction $BAMDIR/$sample"_nodupl_wRG_clip.bam"
mosdepth -t 1 -b 100000 -n $FINALDIR/${sample}_overlap_correction $BAMDIR/$sample"_nodupl_wRG.bam"

### END

### ------------------------------------------------------------------------


### Get average coverage per each window ----------------

cd /scratch/project_2001443/bam/cov_highcovInds_5perSp

#paste per region coverage data for all pixy samples
touch output.txt
touch ind_cov.list
for i in *regions.bed.gz; do
    zcat "$i" | awk -v OFS='\t' '{print $4}' | paste - output.txt > temp_output.txt
    mv temp_output.txt output.txt
    echo "Processed $i"
    echo $i >> ind_cov.list
done

#paste region info as the first column
zcat 108_overlap_correction.regions.bed.gz | awk -v OFS='\t' '{print $1, $2, $3}' > 1stcol.txt #extract region info
paste 1stcol.txt output.txt > mosdepth_highcovinds_100kb_regions_depth.txt
rm 1stcol.txt output.txt

#make the ind_cov.list to have species vcf_id as 1st col and species as 2nd col
nano ind_cov.list
#modify it; resulting in ind_cov_species.list that has 1st row: chr, 2nd row: start, 3rd row: end, followed by ind/species info

#remove regions that are not used: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00 and Scaffold03
grep -v -e mtDNA -e wFhyb -e Spiroplasma -e Scaffold00 -e Scaffold03 mosdepth_highcovinds_100kb_regions_depth.txt > mosdepth_highcovinds_100kb_regions_depth_coreregions.txt ; rm mosdepth_highcovinds_100kb_regions_depth.txt

#get the files  mosdepth_highcovinds_100kb_regions_depth_coreregions.txt, ind_cov_species.list

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/bam/cov_highcovInds_5perSp/mosdepth_highcovinds_100kb_regions_depth_coreregions.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Correlations_whole_genome/depth'
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/bam/cov_highcovInds_5perSp/ind_cov_species.list' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Correlations_whole_genome/depth'







# 4b_bam_coverage_forPIXY.sh

#make a bam list out of all individuals
cd /scratch/project_2001443/barriers_introgr_formica/bam/bam_all
ls *bam > ../bam.list ; cd ..

#extract the samples that are used in PIXY
#Bam list is the PIXY sample list
cd /scratch/project_2001443/barriers_introgr_formica/pixy/depth
INDS=/scratch/project_2001443/barriers_introgr_formica/pixy/groupfiles/group_nonadmixed_minmiss_fbranch_pixy.tab
cut -f1 $INDS > nodupl_pixy_name.list #THIS LIST HAS ONLY SAMPLE NAMES
#wc -l nodupl_pixy_name.list #30 inds

cp /scratch/project_2001443/barriers_introgr_formica/bam/bam_all/bam.list /scratch/project_2001443/barriers_introgr_formica/pixy/depth/bam.pixy.list
nano bam.pixy.list #select manually desired individuals based on the nodupl_pixy_name.list


###
### Compute sequencing depth per sample (after duplicate removal) ------------------------------------------------------------------------
###



#!/bin/bash -l
#SBATCH -J mosdepth_pixy
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/pixy/logs/mosdepth_%j.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/pixy/logs/mosdepth_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-30
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

BAMDIR=/scratch/project_2001443/barriers_introgr_formica/bam/bam_all
FINALDIR=/scratch/project_2001443/barriers_introgr_formica/pixy/depth
FILELIST=/scratch/project_2001443/barriers_introgr_formica/pixy/depth/bam.pixy.list
NAMELIST=/scratch/project_2001443/barriers_introgr_formica/pixy/depth/nodupl_pixy_name_numonly.list

cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILELIST)

# Get sample ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $NAMELIST)

mosdepth -t 1 -b 20000 -n $FINALDIR/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG_clip.bam"
mosdepth -t 1 -b 20000 -n $FINALDIR/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG.bam"

### END

### ------------------------------------------------------------------------


### Get average coverage per each window ----------------

cd /scratch/project_2001443/barriers_introgr_formica/pixy/depth

#paste per region coverage data for all pixy samples
touch output.txt
for i in *regions.bed.gz; do
    zcat "$i" | awk -v OFS='\t' '{print $4}' | paste - output.txt > temp_output.txt
    mv temp_output.txt output.txt
    echo "Processed $i"
done

#paste region info as the first column
zcat 108_overlap_correction.regions.bed.gz | awk -v OFS='\t' '{print $1, $2, $3}' > 1stcol.txt #extract region info
paste 1stcol.txt output.txt > mosdepth_pixy_20kb_regions_depth.txt
rm 1stcol.txt output.txt

#remove regions that are not used: mitchondria (mtDNA), Wolbachia (wFhyb*), and Spiroplasma (Spiroplasma*), as well as Scaffold00 
grep -v -e mtDNA -e wFhyb -e Spiroplasma -e Scaffold00 mosdepth_pixy_20kb_regions_depth.txt > mosdepth_pixy_20kb_regions_depth_coreregions.txt ; rm mosdepth_pixy_20kb_regions_depth.txt


scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/depth/mosdepth_pixy_20kb_regions_depth_coreregions.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy/depth'



################################################### TRIED TO GET THE DEPTH FROM VCF BUT WAS NOT A GOOD IDEA


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
bcftools query -f '[%CHROM:%POS %INFO/DP]\n' $VCF > allsamples_filtered_DP.txt
bcftools query -f '[%CHROM:%POS %INFO/DP]\n' test.vcf > test.txt
bcftools query -f "[%CHROM:%POS %INFO/DP]\n" test.vcf.gz > test.txt
bcftools query -f '[%CHROM:%POS %INFO/DP]' test.vcf.gz > test.txt
bcftools query -f "[%CHROM:%POS %INFO/DP;]" test.vcf.gz | tr ';' '\n' > test.txt

bcftools query -f "[%CHROM:%POS %INFO/DP]\n" --samples "Lai_1w" $VCF > tmp.txt
awk '{gsub(":", " "); print}' tmp.txt > tmp2.txt
bgzip -c tmp2.txt > allsamples_filtered_DP.txt.gz
awk '{gsub(" ", "\t"); print}' tmp2.txt > allsamples_filtered_DP.txt

less -S test.txt


#get these info
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/allsamples_filtered_DP.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy/depth'

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/Formica_hybrid_v1_20kb_regions.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy/depth'





