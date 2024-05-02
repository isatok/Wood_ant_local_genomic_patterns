#!/bin/bash -l
#SBATCH -J twisst_mixedGeo_parentals_aqu_pol_rufa_exs_og_Mi50_sticcs
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_mixedGeo_parentals_aqu_pol_rufa_exs_og_Mi50_sticcs.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_mixedGeo_parentals_aqu_pol_rufa_exs_og_Mi50_sticcs.err
#SBATCH --account=project_2001443
#SBATCH -t 10:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G

cd /scratch/project_2001443/barriers_introgr_formica/twisst

export PATH="/projappl/project_2001443/ete3env/bin:$PATH"

DATADIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs
RESDIR=/scratch/project_2001443/barriers_introgr_formica/twisst/sticcs
GROUPDIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

echo "Running Twisst with sticcs trees"

python3.6 /projappl/project_2001443/twisst/twisst.py \
-t "$DATADIR/phased_exs_out.sticcs.trees.gz" \
-w "$RESDIR/mixedGeo_parentals_exs_og_max50percMi.weights.sticcs.csv.gz" \
-g aqu \
-g rufa \
-g pol \
-g exsecta \
--outgroup exsecta \
--groupsFile $GROUPDIR/group_mixedGeo_parentals_aqu_pol_rufa_sticcs.tab #Fexs is "Fexs_A" in sticcs trees! Inspect why.

###END.
