# 3_RGs.sh

cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl
mkdir /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
ls *bam > nodupl_name.list

#!/bin/bash -l
#SBATCH -J RGs
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/bam/logs/RGs_%j.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/bam/logs/RGs_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-12
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

module load biokit
export PATH="/projappl/project_2001443/picardenv/bin:$PATH"  #use picard v 2.21.4, which was used for my other data as well 

FINALDIR=/scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG

cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p nodupl_name.list)

# Get sample ID
sample=${file%_nodupl*}


###
### Add read group
###

# Add RGs
picard AddOrReplaceReadGroups \
    I=$file \
    O=$FINALDIR/$sample"_nodupl_wRG.bam" \
    RGID=$sample \
    RGPL=illumina \
    RGLB=1 \
    RGSM=$sample \
    RGPU=1 \
    TMP_DIR=/scratch/project_2001443/tmp

# Index
samtools index $FINALDIR/$sample"_nodupl_wRG.bam"


### END
