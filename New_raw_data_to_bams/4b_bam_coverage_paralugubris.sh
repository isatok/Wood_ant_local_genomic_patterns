# 4b_bam_coverage_paralugubris.sh

###
### Compute sequencing depth per sample (after duplicate removal) ------------------------------------------------------------------------
###


cd /scratch/project_2001443/paralugubris/bam/nodupl
ls *bam > ../bam.list ; cd ..


#!/bin/bash -l
#SBATCH -J mosdepth_paralugubris
#SBATCH -o /scratch/project_2001443/paralugubris/bam/logs/mosdepth_%j.out
#SBATCH -e /scratch/project_2001443/paralugubris/bam/logs/mosdepth_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 02:00:00
#SBATCH -p small
#SBATCH --array=1-2
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

FINALDIR=/scratch/project_2001443/paralugubris/bam/nodupl_RG

cd /scratch/project_2001443/paralugubris/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p nodupl_name.list)

# Get sample ID
sample=${file%_nodupl*}

mosdepth -t 1 -b 10000 -n -x ../stats/coverage/$sample $FINALDIR/$sample"_nodupl_wRG.bam"
mosdepth -t 1 -b 10000 -n ../stats/coverage/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG.bam"

### END


### Get coverage for overlap corrected data  --------------------

cd /scratch/project_2001443/paralugubris/bam/stats/coverage
rm -rf *global*

for i in *overlap_correction.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | awk '{sum+=$4;} END{print sum/NR;}'
done > paralugubris.overlap_correction.mosdepth.summary.tmp

grep -v 'txt' paralugubris.overlap_correction.mosdepth.summary.tmp > tmp1
grep 'txt' paralugubris.overlap_correction.mosdepth.summary.tmp | perl -npe 's/_overlap_correction.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > paralugubris.overlap_correction.mosdepth.summary.txt ; rm *tmp*

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/paralugubris/bam/stats/coverage/paralugubris.overlap_correction.mosdepth.summary.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/mapping_and_insert_stats/'

### Get coverage for *NON* overlap corrected data  --------------------

for i in s???.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | awk '{sum+=$4;} END{print sum/NR;}'
done > paralugubris.no_overlap_correction.mosdepth.summary.tmp

grep -v 'txt' paralugubris.no_overlap_correction.mosdepth.summary.tmp > tmp1
grep 'txt' paralugubris.no_overlap_correction.mosdepth.summary.tmp | perl -npe 's/.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > paralugubris.no_overlap_correction.mosdepth.summary.txt ; rm *tmp*

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/paralugubris/bam/stats/coverage/paralugubris.no_overlap_correction.mosdepth.summary.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/mapping_and_insert_stats/'

scp puhti:/scratch/project_2001443/paralugubris/bam/coverage/paralugubris.overlap_correction.mosdepth.summary.txt .

# in R
# options(stringsAsFactors=F)
# tt = read.table("paralugubris.overlap_correction.mosdepth.summary.txt", h = F)
# summary(tt$V2)
# overall mean: XXX (XX - XX)
