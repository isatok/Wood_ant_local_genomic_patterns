# 4_bam_coverage.sh

###
### Compute sequencing depth per sample (after duplicate removal) ------------------------------------------------------------------------
###


cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl
ls *bam > ../bam.list ; cd ..


#!/bin/bash -l
#SBATCH -J mosdepth
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/bam/logs/mosdepth_%j.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/bam/logs/mosdepth_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-12
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

FINALDIR=/scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG

cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p nodupl_name.list)

# Get sample ID
sample=${file%_nodupl*}

mosdepth -t 1 -b 10000 -n -x ../stats/coverage/$sample $FINALDIR/$sample"_nodupl_wRG.bam"
mosdepth -t 1 -b 10000 -n ../stats/coverage/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG.bam"

### END


### Get coverage for overlap corrected data  --------------------

cd /scratch/project_2001443/barriers_introgr_formica/bam/stats/coverage
rm -rf *global*

for i in *overlap_correction.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | awk '{sum+=$4;} END{print sum/NR;}'
done > all.overlap_correction.mosdepth.summary.tmp

grep -v 'txt' all.overlap_correction.mosdepth.summary.tmp > tmp1
grep 'txt' all.overlap_correction.mosdepth.summary.tmp | perl -npe 's/_overlap_correction.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > all.overlap_correction.mosdepth.summary.txt ; rm *tmp*

#this makes the right computation (excludes scaffs 00 & 03), although it doesn't take into account the different chr lengths which is åperhaps now ok
for i in *overlap_correction.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | grep -v 'region' | grep -v 'Scaffold03' | grep -v 'Scaffold00' | awk '{sum+=$4;} END{print sum/NR;}'
done > all.overlap_correction_no0003.mosdepth.summary.tmp
grep -v 'txt' all.overlap_correction_no0003.mosdepth.summary.tmp > tmp1
grep 'txt' all.overlap_correction_no0003.mosdepth.summary.tmp | perl -npe 's/_overlap_correction_no0003.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > all.overlap_correction_no0003.mosdepth.summary.txt ; rm *tmp*

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/bam/bam_all/stats/coverage/all.overlap_correction_no0003.mosdepth.summary.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/mapping_and_insert_stats/'

### Get coverage for *NON* overlap corrected data  --------------------

for i in RN???.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | awk '{sum+=$4;} END{print sum/NR;}'
done > all.no_overlap_correction.mosdepth.summary.tmp

grep -v 'txt' all.no_overlap_correction.mosdepth.summary.tmp > tmp1
grep 'txt' all.no_overlap_correction.mosdepth.summary.tmp | perl -npe 's/.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > all.no_overlap_correction.mosdepth.summary.txt ; rm *tmp*

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/bam/stats/coverage/all.no_overlap_correction.mosdepth.summary.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/mapping_and_insert_stats/'

scp puhti:/scratch/project_2001443/barriers_introgr_formica/bam/coverage/all.overlap_correction.mosdepth.summary.txt .

### MODIFY THIS WHEN DECIDING WHICH SAMPLES TO KEEP (ESP #217 MEAN COV=2.16 IS LOW ###
# in R
# options(stringsAsFactors=F)
# tt = read.table("all.overlap_correction.mosdepth.summary.txt", h = F)
# summary(tt$V2)
# overall mean: XXX (XX - XX)
