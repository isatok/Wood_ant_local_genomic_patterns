###
### Compute sequencing depth per sample (after duplicate removal) ------------------------------------------------------------------------
###
WORK WITH THIS!!


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

mosdepth -t 1 -b 10000 -n -x ../stats/coverage/$sample $FINALDIR/$sample"_nodupl_wRG.bam"
mosdepth -t 1 -b 10000 -n ../stats/coverage/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG.bam"

