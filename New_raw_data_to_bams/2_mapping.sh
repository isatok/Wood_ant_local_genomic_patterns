#!/bin/bash -l
#SBATCH -J map
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/bam/logs/map_%j.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/bam/logs/map_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-12
#SBATCH --ntasks 8
#SBATCH --mem=12G
#SBATCH --mail-user=ina.satokangas@helsinki.fi
#SBATCH --mail-type=END


###
### 0. Prep -------------------------------------------------------------------
###

# load modules
module load biokit
module load r-env
export PATH="/projappl/project_2001443/picardenv/bin:$PATH"  #use picard v 2.21.4, which was used for my other data as well 

# Set work directory
cd /scratch/project_2001443/barriers_introgr_formica/

# Define paths
REFPATH=/scratch/project_2001443/reference_genome
FASTQPATH=/scratch/project_2001443/barriers_introgr_formica/fastq/trim
BAMPATH=/scratch/project_2001443/barriers_introgr_formica/bam

# Get file name and sample ID
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p fastq/name.list)
shortfile=${file}

echo "###### STARTING file: $file"


###
### 1. Map, sort & index ------------------------------------------------------
###

echo "###### MAPPING"

bwa mem -t8 $REFPATH/Formica_hybrid_v1_wFhyb_Sapis.fa $FASTQPATH/$file"_1.trim.pair.fq.gz" $FASTQPATH/$file"_2.trim.pair.fq.gz" | samtools sort -@8 -m512M -o $BAMPATH/raw/${shortfile}".bam" -

samtools index -@4 $BAMPATH/raw/${shortfile}".bam"


###
### 2. Compute insert size distribution ---------------------------------------
###

echo "###### COLLECTING INSERT SIZES"

picard CollectInsertSizeMetrics \
I=$BAMPATH/raw/${shortfile}".bam" \
O=$BAMPATH/stats/${shortfile}"_insert_size_metrics.txt" \
H=$BAMPATH/stats/${shortfile}"_insert_size_hist.pdf"


###
### 3. Filter duplicates ----------------------------------------------------
###

echo "###### FILTERING DUPLICATES"

picard MarkDuplicates \
I=$BAMPATH/raw/${shortfile}".bam" \
O=$BAMPATH/nodupl/${shortfile}"_nodupl.bam" \
M=$BAMPATH/nodupl/stats/${shortfile}"_dupl_metrics.txt" \
REMOVE_DUPLICATES=T


###
### 4. Index and compute stats ------------------------------------------------
###

echo "###### INDEX & GET FLAGSTATS"

samtools index -@4 $BAMPATH/nodupl/${shortfile}"_nodupl.bam" && \
samtools flagstat -@4 $BAMPATH/nodupl/${shortfile}"_nodupl.bam" > $BAMPATH/nodupl/stats/${shortfile}"_nodupl.flagstat" && \
samtools idxstats -@4 $BAMPATH/nodupl/${shortfile}"_nodupl.bam" > $BAMPATH/nodupl/stats/${shortfile}"_nodupl.idxstat"

echo "###### DONE! sample ID: $file"

### This is the end.
