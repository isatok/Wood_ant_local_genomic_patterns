#update Simon's scripts if needed: 
cd /scratch/project_2001443/analysis/genomics_simon
git clone https://github.com/simonhmartin/genomics_general.git

---

cd /scratch/project_2001443/vcf/geno

sinteractive --account project_2001443 --mem 6000

DATADIR=/scratch/project_2001443/vcf/geno
RESDIR=/scratch/project_2001443/analysis/twisst/trees/exs
GENO=phased.exsecta.geno.gz
PLOIDYTAB=/scratch/project_2001443/vcf/ind_lists/ploidy.tab

export PYTHONPATH="/scratch/project_2001443/analysis/genomics_simon/genomics_general:$PYTHONPATH"
export PATH="/projappl/project_2001443/ete3env/bin:$PATH" #This contains NUMPY

  
  python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/filterGenotypes.py \
  -t 4 \
  -i $DATADIR/$GENO \
  --ploidyFile /scratch/project_2001443/vcf/ind_lists/ploidy.tab \
  --outputGenoFormat bases \
  -o $DATADIR/phased.exsecta.bases.geno.gz
 
 python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/freq.py \
 -g $DATADIR/phased.exsecta.bases.geno.gz --target derived --threads 10 --asCounts TRUE \
 -o $DATADIR/phased.exsecta.bases.csv --popsFile

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

#>>>>>>>>>>>>>>>>>>>>>>INA CONTINUE  FROM HERE 16.08.2023<<<<<<<<<<<<<<<<<<<<<<<<<<<

sinteractive --account project_2001443 --mem 6000

cd /scratch/project_2001443/analysis/twisst/trees_simon/

INPUT=/scratch/project_2001443/vcf/geno/phased.exsecta.bases.csv
REFERENCE=/scratch/project_2001443/reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa.fai

python3 /scratch/project_2001443/analysis/twisst/trees_simon/sticcs.py \
-i $INPUT -o test_160823 -l $REFERENCE --allowSecondChances
