# 5_clip_overlaps.sh


# Create the input file list
cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
ls *bam | sort -nk1,1 > input.list ; cd .. #sort numeric, according to position 1,1

# Create the output folder
mkdir /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip


###
### sbatch script below ------------------------------------------------------------------------
###

####MODIFY AND RUN FROM HERE ONWARDS ###

#!/bin/bash -l
#SBATCH -J clipOverlapRG
#SBATCH -o /scratch/project_2001443/bam/logs/clipOverlap_RG_%j.out
#SBATCH -e /scratch/project_2001443/bam/logs/clipOverlap_RG_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-71
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

# Load modules
module load biokit
module load bioconda/3
source activate my_seqdata

# Define directories
cd /scratch/project_2001443/bam/nodupl_RG
FINALDIR=/scratch/project_2001443/bam/nodupl_RG_clip

# Get file & sample ID
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p input.list)
sample=${file%_wRG*}

echo "### Processing $file"


###
### Clip overlaps
###

# Reads are sorted by coordinates (see "samtools sort" command after mapping)
bam clipOverlap --in $FINALDIR/${sample}"_wRG.bam" --out $FINALDIR/${sample}"_wRG_clip.bam" --stats --params


###
### Index
###

samtools index $FINALDIR/$sample"_wRG_clip.bam"



### END
