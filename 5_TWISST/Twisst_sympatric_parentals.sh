### TWISST is not used for the current manuscript (updated in May 2024) ###

#!/bin/bash -l
#SBATCH -J twisst_symparsamples_exsog_Mi50
#SBATCH --account=project_2001443
#SBATCH -t 10:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G

cd /scratch/project_2001443/analysis/twisst
export PATH="/projappl/project_2001443/ete3env/bin:$PATH"


for x in 50; do

echo "Running Twisst for window size $x"

python3.6 /projappl/project_2001443/twisst/twisst.py \
-t "trees/exs/phased_exs_out_max50percMi.phyml_bionj.w"$x".trees.gz" \
-w "results/parentals/exs_outgroup/symparsamples_exs_og_max50percMi.w"$x".weights.csv.gz" \
-g pratensis_fi \
-g aquilonia_fi \
-g rufa_fi \
-g lugubris_fi \
-g exsecta \
--outgroup exsecta \
--groupsFile group_symparsamples.tab

done
twisst_symparsamples_exsog_50Mi.sh (END)
