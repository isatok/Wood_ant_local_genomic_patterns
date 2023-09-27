
###
### 1. Prepare group files for PIXY
###

# see at the bottom (takes a lot of space so instructions to make these are at the bottom of the script) #

###
### 2. Compute PIXY
###

#!/bin/bash -l
#SBATCH -J pixy_langholmen
#SBATCH -o /scratch/project_2001443/analysis/pixy_langholmen/logs/pixy1.out
#SBATCH -e /scratch/project_2001443/analysis/pixy_langholmen/logs/pixy1.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

cd /scratch/project_2001443/analysis/pixy_langholmen/

pixy --stats pi \
--vcf /scratch/project_2001443/analysis/pixy_langholmen/vcf/allnoDiploidMales_filtered.vcf.gz \
--populations pop.list \
--window_size 10000 \
--output_prefix allnoDiploidMales \
--n_cores 8


# sort output - FST
head -n1 allFems_fst.txt > tmp
grep -v 'pop' allFems_fst.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_fst_sort.txt

# sort output - Dxy
head -n1 allFems_dxy.txt > tmp
grep -v 'pop' allFems_dxy.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_dxy_sort.txt

# sort output - Pi
head -n1 allFems_pi.txt > tmp
grep -v 'pop' allFems_pi.txt | sort -k2,2 -k3,3n > tmp2
cat tmp tmp2 > allFems_pi_sort.txt

rm tmp tmp2

###This is the end




###
### **. Prepare group files for PIXY
###

# see at the bottom (takes a lot of space so instructions to make these are at the bottom of the script) #
FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab

cut -f3 $FULLSAMPLE | sort | uniq   #check which groups are possible (w-out geogr. information)
# adm1
# adm2
# aqu
# exsecta
# lug
# paral
# pol
# prat
# rufa

cut -f4 $FULLSAMPLE | sort | uniq   #check which groups are possible (WITH geogr. information)
# adm1_ceu
# adm1_fi
# adm2_ceu
# adm2_fi
# aqu_ceu
# aqu_fi
# exsecta
# lug_fi
# paral_ceu
# pol_ceu
# prat_fi
# rufa_ceu
# rufa_fi

#Borrow the group files created for global Fst estimation; for record how they are created:
cd /scratch/project_2001443/barriers_introgr_formica/fst_global/groupfiles

# F. rufa
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa" {print $1}' > rufa_all.tab #12 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_fi" {print $1}' > rufa_fi.tab #6 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_ceu" {print $1}' > rufa_ceu.tab #6 inds
shuf -n 3 rufa_fi.tab > rufa_all_6inds.tab; shuf -n 3 rufa_ceu.tab >> rufa_all_6inds.tab #6 random inds (3fi, 3 ceu)
shuf -n 3 rufa_fi.tab > rufa_all_6inds_iter2.tab; shuf -n 3 rufa_ceu.tab >> rufa_all_6inds_iter2.tab #6 random inds (3fi, 3 ceu), 2nd iteration

# F. aquilonia
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" {print $1}' > aquilonia_all.tab #45 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_fi" {print $1}' > aquilonia_fi.tab #39 inds
shuf -n 6 aquilonia_fi.tab > aquilonia_fi_6inds.tab #6 random inds
shuf -n 6 aquilonia_fi.tab > aquilonia_fi_6inds_iter2.tab #6 random inds, 2nd iteration
shuf -n 3 aquilonia_fi.tab > aquilonia_fi_3inds.tab #3 random inds
shuf -n 3 aquilonia_fi.tab > aquilonia_fi_3inds_iter2.tab #3 random inds, 2nd iteration
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_ceu" {print $1}' > aquilonia_ceu.tab #6 inds
shuf -n 3 aquilonia_fi.tab > aquilonia_all_6inds.tab; shuf -n 3 aquilonia_ceu.tab >> aquilonia_all_6inds.tab #6 inds
shuf -n 3 aquilonia_fi.tab > aquilonia_all_6inds_iter2.tab; shuf -n 3 aquilonia_ceu.tab >> aquilonia_all_6inds_iter2.tab #6 random inds, 2nd iteration

cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "Lai_1w" || $1 == "Lai_2w" || $1 == "Loa_1w" {print $1}' > aquilonia_scot.tab #3 inds
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CBAQ1_1w" || $1 == "CBAQ2_2w" || $1 == "CBAQ3_1w" {print $1}' > aquilonia_swi.tab #3 inds

# F. polyctena
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $1}' > polyctena_ceu.tab #6 inds

GROUPFILE=group_parentals_aqu_pol_rufa_pixy.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" || \
$2 == "rufa" || $2 == "pol" || $2 == "exsecta" {print $0}' > $GROUPFILE.tmp


