#!/bin/bash -l
#SBATCH -J twisst_mixedGeo_parentals_aqu_pol_rufa_exs_og_Mi50_phyml
#SBATCH --account=project_2001443
#SBATCH -t 10:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G

#####CONTINUE EDITING FROM HERE FOR 2023 AUTUMN ANALYSES #####


cd /scratch/project_2001443/analysis/twisst

module load bioconda/3
source activate ete3


for x in 50; do

echo "Running Twisst for window size $x"

python3.6 /projappl/project_2001443/twisst/twisst.py \
-t "trees/exs/phased_exs_out_max50percMi.phyml_bionj.w"$x".trees.gz" \
-w "results/parentals/exs_outgroup/parentals_fi_exs_og_max50percMi.w"$x".weights.csv.gz" \
-g pratensis_fi \
-g aquilonia_fi \
-g rufa_fi \
-g lugubris_fi \
-g exsecta \
--outgroup exsecta \
--groupsFile group_parentals_fi.tab

done

twisst_parentals_fi_exsog_50Mi.sh (END)
