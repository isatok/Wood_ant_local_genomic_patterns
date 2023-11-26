### Flagstats

#!/bin/bash -l
#SBATCH -J flagstats
#SBATCH -o /scratch/project_2001443/bam/logs/flags_%j.out
#SBATCH -e /scratch/project_2001443/bam/logs/flags_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 06:00:00
#SBATCH -p small
#SBATCH --array=2-20
#SBATCH --ntasks 8
#SBATCH --mem=12G

module load biokit

cd /scratch/project_2001443
BAMPATH=/scratch/project_2001443/bam

file=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2001443/bam/nodupl_RG_clip/pierre.list)

samtools flagstat -@4 $BAMPATH/nodupl_RG_clip/$file > $BAMPATH/stats/map_stats/$file"_nodupl.flagstat"
/scratch/project_2001443/bam/pierre_flagstat.sh

### Add read groups and compute coverage

#!/bin/bash -l
#SBATCH -J RG_coverage
#SBATCH -o /scratch/project_2001443/bam/logs/RG_coverage_%j.out
#SBATCH -e /scratch/project_2001443/bam/logs/RG_coverage_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-71
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

module load biokit

FINALDIR=/scratch/project_2001443/bam/nodupl_RG

cd /scratch/project_2001443/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p nodupl_name.list)

# Get sample ID
sample=${file%_nodupl*}


###
### Add read group
###

# Add RGs
java -Xmx4G -jar /appl/soft/bio/picard/picard-tools-2.21.4/picard.jar AddOrReplaceReadGroups \
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


###
### Compute coverage per sample
###

module load bioconda/3
source activate myenv

mosdepth -t 1 -b 10000 -n -x ../stats/coverage/$sample $FINALDIR/$sample"_nodupl_wRG.bam"
mosdepth -t 1 -b 10000 -n ../stats/coverage/${sample}_overlap_correction $FINALDIR/$sample"_nodupl_wRG.bam"

###
### Compute coverage for Scaffold03 only
###

#!/bin/bash -l
#SBATCH -J coverage_scaff3
#SBATCH -o /scratch/project_2001443/bam/logs/coverage_scaff3_%j.out
#SBATCH -e /scratch/project_2001443/bam/logs/coverage_scaff3_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-91
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

module load biokit
BAMDIR=/scratch/project_2001443/bam/nodupl_RG_clip

cd /scratch/project_2001443/bam/nodupl_RG_clip

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p bam.list)

# Get sample ID
sample=${file%_nodupl*}

###
### Compute coverage per sample
###

module load bioconda/3
source activate my_seqdata

mosdepth -t 1 -b 100000 -c Scaffold03  -n -x ../stats/coverage_scaff03/$sample $BAMDIR/$sample"_nodupl_wRG_clip.bam"
mosdepth -t 1 -b 100000 -c Scaffold03 -n ../stats/coverage_scaff03/${sample}_overlap_correction $BAMDIR/$sample"_nodupl_wRG_clip.bam"
