
### Not used for the current manuscript (updated in May 2024) ###

#update Simon's scripts if needed: 
cd /scratch/project_2001443/analysis/genomics_simon
git clone https://github.com/simonhmartin/genomics_general.git

---

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

sinteractive --account project_2001443 --mem 6000

DATADIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
GENO=phased.exsecta.geno.gz
PLOIDYTAB=ploidy.tab

export PYTHONPATH="/scratch/project_2001443/analysis/genomics_simon/genomics_general:$PYTHONPATH"
#export PATH="/projappl/project_2001443/ete3env/bin:$PATH" #This contains NUMPY
export PATH="/projappl/project_2001443/phymlenv/bin:$PATH"

### Transform the data from the geno.gz format into bases.geno.gz, which has each sample divided into two columns A&B corresponding to the phased haplotypes
  
  python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/filterGenotypes.py \
  -t 4 \
  -i $DATADIR/$GENO \
  --ploidyFile $DATADIR/$PLOIDYTAB \
  --outputGenoFormat bases \
  -o $DATADIR/phased.exsecta.bases.geno.gz

### Transform the data from the bases.geno.gz format into bases.csv, which (in this case) gives information on whether each allele in each individual is derived or ancestral
### -> this results with 0 and 1 values (corresponding to ancestral and derived with respect to the exsecta sample).
 
python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/freq.py \
-g $DATADIR/phased.exsecta.bases.geno.gz -o $DATADIR/phased.exsecta.bases.csv \
--indFreqs --target derived  --threads 10  --asCounts --verbose --ploidyFile ploidyAB.tab


 
#### build trees with sticcs.py

#input file: /scratch/project_2001443/vcf/geno/phased.exsecta.bases.csv
#needs imputed or phased data (since it allows no missing data)
#sensitive to genotyping (and imputation) errors. when zooming in to narrow regions this manifests as spurious trees but with broader distances works well
#Simon wants to see:
#-relative proportions of topologies genome-wide, and patterns along the chromosomes, in comparison to PhyML trees
#-also interesting to make a comparison with relate & tsinfer but not priority now
#download the script from discord

#python3 /scratch/project_2001443/analysis/twisst/trees_simon/sticcs.py \
#-i input.freqs.tsv.gz -o output_prefix -l reference.fai
#one more option is --allowSecondChances

sinteractive --account project_2001443 --mem 6000

#export PATH="/projappl/project_2001443/phymlenv/bin:$PATH" # No need to export if already did it during the same session (above)

cd /scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs

DATADIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
RESDIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs
BASES=phased.exsecta.bases.csv
REFERENCE=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa.fai
OUTFILE=phased_exs_out.sticcs

python3 /scratch/project_2001443/analysis/twisst/trees_simon/sticcs.py \
-i $DATADIR/$BASES -o $RESDIR/$OUTFILE -l $REFERENCE --allowSecondChances


#####------------------- sticcs as a bash script

#!/bin/bash -l
#SBATCH -J sticcs_trees_exs
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs/logs/sticcs_trees_exs.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs/logs/sticcs_trees_exs.err
#SBATCH --account=project_2001443
#SBATCH -t 08:00:00
#SBATCH -p small
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-type=END

cd /scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs

export PATH="/projappl/project_2001443/phymlenv/bin:$PATH"

DATADIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
RESDIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs
BASES=phased.exsecta.bases.csv
REFERENCE=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa.fai
OUTFILE=phased_exs_out.sticcs

python3 /scratch/project_2001443/analysis/twisst/trees_simon/sticcs.py \
-i $DATADIR/$BASES -o $RESDIR/$OUTFILE -l $REFERENCE --allowSecondChances


###END.
