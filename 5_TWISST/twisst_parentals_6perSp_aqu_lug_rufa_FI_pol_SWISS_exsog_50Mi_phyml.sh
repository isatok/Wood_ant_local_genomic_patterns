#!/bin/bash -l
#SBATCH -J twisst_parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS_exsog_50Mi
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS_exsog_50Mi_phyml.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/twisst/logs/twisst_parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS_exsog_50Mi_phyml.err
#SBATCH --account=project_2001443
#SBATCH -t 5:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --mail-type=END

cd /scratch/project_2001443/barriers_introgr_formica/twisst/
export PATH="/projappl/project_2001443/ete3env/bin:$PATH"

DATADIR=/scratch/project_2001443/barriers_introgr_formica/local_trees/phyml
RESDIR=/scratch/project_2001443/barriers_introgr_formica/twisst/phyml
GROUPDIR=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

for x in 50; do

echo "Running Twisst for window size $x"

python3.6 /projappl/project_2001443/twisst/twisst.py \
-t "$DATADIR/phased_exs_out_max50percMi.phyml_bionj.w"$x".trees.gz" \
-w "$RESDIR/parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS_exs_og_max50percMi.w"$x".weights.csv.gz" \
-g aquilonia \
-g rufa \
-g lugubris_fi \
-g polyctena_ceu \
-g exsecta \
--outgroup exsecta \
--groupsFile $GROUPDIR/group_parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS.tab
done

# twisst_parentals_6perSp_aqu_lug_rufa_FI_pol_SWISS_exsog_50Mi_phyml.sh (END)
