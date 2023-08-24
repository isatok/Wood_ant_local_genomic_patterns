# 5_clip_overlaps.sh


# Create the input file list
cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
ls *bam | sort -nk1,1 > input.list ; cd .. #sort numeric, according to position 1,1

# Create the output folder
mkdir /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip


###
### sbatch script below ------------------------------------------------------------------------
###

#!/bin/bash -l
#SBATCH -J clipOverlap
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/bam/logs/clipOverlap_%j.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/bam/logs/clipOverlap_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-12
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

# Load modules
module load biokit
export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

# Define directories
cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
INPUTDIR=/scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
FINALDIR=/scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip

# Get file & sample ID
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p input.list)
sample=${file%_nodupl_wRG*}

echo "### Processing $file"


###
### Clip overlaps
###

# Reads are sorted by coordinates (see "samtools sort" command after mapping)
bam clipOverlap --in $INPUTDIR/${sample}"_nodupl_wRG.bam" --out $FINALDIR/${sample}"_nodupl_wRG_clip.bam" --stats --params


###
### Index
###

samtools index $FINALDIR/$sample"_nodupl_wRG_clip.bam"



### END

### Check that the bam files are alright: -------
#(Despite the errors, last project (Satokangas et al 2023) worked well, so no reason for now to freak out about these.)

## Unclipped ## -------
/scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG
picard ValidateSamFile I=RN415_nodupl_wRG.bam R=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa MODE=SUMMARY

## HISTOGRAM	java.lang.String
#Error Type	Count
#ERROR:INVALID_TAG_NM	26664

## Clipped ## -------
cd /scratch/project_2001443/barriers_introgr_formica/bam/nodupl_RG_clip
picard ValidateSamFile I=RN415_nodupl_wRG_clip.bam R=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa MODE=SUMMARY

## HISTOGRAM	java.lang.String
#Error Type	Count
#ERROR:CIGAR_MAPS_OFF_REFERENCE	655
#ERROR:INVALID_CIGAR	1055424
#ERROR:INVALID_TAG_NM	5704471
#ERROR:MISMATCH_MATE_CIGAR_STRING	12683598
