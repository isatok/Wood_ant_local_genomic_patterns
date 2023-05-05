cd /scratch/project_2001443/vcf/geno

sinteractive ...

DATADIR=/scratch/project_2001443/vcf/geno
RESDIR=/scratch/project_2001443/analysis/twisst/trees/exs
GENO=phased.exsecta.geno.gz
PLOIDYTAB=/scratch/project_2001443/vcf/ind_lists/ploidy.tab

export PYTHONPATH="/scratch/project_2001443/analysis/genomics_simon/genomics_general:$PYTHONPATH"
export PATH="/projappl/project_2001443/ete3env/bin:$PATH" #This includes NUMPY

  
  python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/filterGenotypes.py \
  -t 4 \
  -i $DATADIR/$GENO \
  --ploidyFile /scratch/project_2001443/vcf/ind_lists/ploidy.tab \
  --outputGenoFormat bases \
  -o $DATADIR/phased.exsecta.bases.geno.gz
   
 ## next run freq.py for this  $DATADIR/phased.exsecta.bases.geno.gz file. find out about how to do it.
 
 python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/freq.py \
 -g $DATADIR/phased.exsecta.bases.geno.gz --target derived --threads 10 --asCounts TRUE \
 -o $DATADIR/phased.exsecta.bases.csv --popsFile

#MODIFY EG PLOIDYTAB SO IT IS AS POPSFILE (SAMPLENAME SAMPLENAME)
