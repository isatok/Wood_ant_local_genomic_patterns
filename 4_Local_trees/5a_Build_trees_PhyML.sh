#Build local trees with PhyML, allowing max 50% of missing data per individual, for 50 SNPs per window.
#Use F. exsecta as the outgroup, include all sequenced samples excluding 110 (collaborative sample not published here) and 105 (too much missing data).

#!/bin/bash -l
#SBATCH -J phyml_trees_exs_Mi50
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/local_trees/phyml/logs/phyml_trees_exs_Mi50.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/local_trees/phyml/logs/phyml_trees_exs_Mi50.err
#SBATCH --account=project_2001443
#SBATCH -t 08:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --mail-type=END

cd /scratch/project_2001443/barriers_introgr_formica/local_trees/phyml

DATADIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
RESDIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/phyml

export PATH="/projappl/project_2001443/phymlenv/bin:$PATH" 
export PYTHONPATH="/scratch/project_2001443/analysis/genomics_simon/genomics_general:$PYTHONPATH"

# Fraction of missing data allowed
fractMiss=50

for x in 50 ; do

# compute max missing data per individual allowed (set to 50%)
(( Mi = x - x*fractMiss/100 ))

echo "Inferring trees with window size $x SNPs, at least $Mi sites per individual"

python3 /scratch/project_2001443/analysis/genomics_simon/genomics_general/phylo/phyml_sliding_windows.py -T 4 \
  -g $DATADIR/phased.exsecta.geno.gz \
  --prefix $RESDIR/phased_exs_out_max${fractMiss}percMi.phyml_bionj.w$x \
  -w $x -Mi $Mi --windType sites --model GTR --optimise n

done

###END.

trees_exs_Mi50.sh (END)
