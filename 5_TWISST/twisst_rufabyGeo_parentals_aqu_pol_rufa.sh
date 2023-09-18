#!/bin/bash -l
#SBATCH -J twisst_rufabyGeo_parentals_aqu_pol_rufa_exs_og_Mi50_phyml
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_rufabyGeo_parentals_aqu_pol_rufa_exs_og_Mi50_phyml.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_rufabyGeo_parentals_aqu_pol_rufa_exs_og_Mi50_phyml.err
#SBATCH --account=project_2001443
#SBATCH -t 10:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G

cd /scratch/project_2001443/barriers_introgr_formica/twisst

export PATH="/projappl/project_2001443/ete3env/bin:$PATH"

DATADIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/phyml
RESDIR=/scratch/project_2001443/barriers_introgr_formica/twisst/phyml
GROUPDIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

for x in 50; do

echo "Running Twisst for window size $x"

python3.6 /projappl/project_2001443/twisst/twisst.py \
-t "$DATADIR/phased_exs_out_max50percMi.phyml_bionj.w"$x".trees.gz" \
-w "$RESDIR/rufabyGeo_parentals_exs_og_max50percMi.w"$x".weights.csv.gz" \
-g aqu \
-g rufa_fi \
-g rufa_ceu \
-g pol \
-g exsecta \
--outgroup exsecta \
--groupsFile $GROUPDIR/group_rufabyGeo_parentals_aqu_pol_rufa.tab
done

twisst_rufabyGeo_parentals_aqu_pol_rufa.sh (END)
