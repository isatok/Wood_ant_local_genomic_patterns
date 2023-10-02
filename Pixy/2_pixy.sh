#### This script is to compute local Pi, Dxy and Fst estimates with PIXY ####

mkdir /scratch/project_2001443/barriers_introgr_formica/pixy
mkdir /scratch/project_2001443/barriers_introgr_formica/pixy/groupfiles
mkdir /scratch/project_2001443/barriers_introgr_formica/pixy/logs

cd /scratch/project_2001443/barriers_introgr_formica/pixy


###
### 1. Prepare group files for PIXY ------
###

# see at the bottom (takes a lot of space so instructions to make these are at the bottom of the script) #


###
### 2. Compute PIXY ------
###

#!/bin/bash -l
#SBATCH -J pixy
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/pixy/logs/pixy.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/pixy/logs/pixy.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

cd /scratch/project_2001443/barriers_introgr_formica/pixy/

pixy --stats pi dxy fst \
--vcf /scratch/project_2001443/barriers_introgr_formica/gvfc/allsamples_filtered.vcf.gz \
--populations groupfiles/group_parentals_all_pixy.tab \
--window_size 20000 \
--output_prefix allspecies_balanced \
--n_cores 8


# sort output - FST
head -n1 allspecies_balanced_fst.txt > tmp
grep -v 'pop' allspecies_balanced_fst.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allspecies_balanced_fst_sort.txt

# sort output - Dxy
head -n1 allspecies_balanced_dxy.txt > tmp
grep -v 'pop' allspecies_balanced_dxy.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allspecies_balanced_dxy_sort.txt

# sort output - Pi
head -n1 allspecies_balanced_pi.txt > tmp
grep -v 'pop' allspecies_balanced_pi.txt | sort -k2,2 -k3,3n > tmp2
cat tmp tmp2 > allspecies_balanced_pi_sort.txt

rm tmp tmp2

###This is the end


#get the results for plotting
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/pixy/*.txt' '/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Pixy'


###
### **. Prepare group files for PIXY ------
###

### Check the possible groups

FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab

# cut -f3 $FULLSAMPLE | sort | uniq   #check which groups are possible (w-out geogr. information)
# # adm1
# # adm2
# # aqu
# # exsecta
# # lug
# # paral
# # pol
# # prat
# # rufa

# cut -f4 $FULLSAMPLE | sort | uniq   #check which groups are possible (WITH geogr. information)
# # adm1_ceu
# # adm1_fi
# # adm2_ceu
# # adm2_fi
# # aqu_ceu
# # aqu_fi
# # exsecta
# # lug_fi
# # paral_ceu
# # pol_ceu
# # prat_fi
# # rufa_ceu
# # rufa_fi


### Borrow the group files created for global Fst estimation; for record how they are created:

# cd /scratch/project_2001443/barriers_introgr_formica/fst_global/groupfiles
# # F. rufa
# cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa" {print $1}' > rufa_all.tab #12 inds
# cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_fi" {print $1}' > rufa_fi.tab #6 inds
# cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_ceu" {print $1}' > rufa_ceu.tab #6 inds
# shuf -n 3 rufa_fi.tab > rufa_all_6inds.tab; shuf -n 3 rufa_ceu.tab >> rufa_all_6inds.tab #6 random inds (3fi, 3 ceu)   <<< THIS IS USED FOR THE FIRST TRY <<<
# shuf -n 3 rufa_fi.tab > rufa_all_6inds_iter2.tab; shuf -n 3 rufa_ceu.tab >> rufa_all_6inds_iter2.tab #6 random inds (3fi, 3 ceu), 2nd iteration

# # F. aquilonia
# cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" {print $1}' > aquilonia_all.tab #45 inds
# cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_fi" {print $1}' > aquilonia_fi.tab #39 inds
# shuf -n 6 aquilonia_fi.tab > aquilonia_fi_6inds.tab #6 random inds
# shuf -n 6 aquilonia_fi.tab > aquilonia_fi_6inds_iter2.tab #6 random inds, 2nd iteration
# shuf -n 3 aquilonia_fi.tab > aquilonia_fi_3inds.tab #3 random inds
# shuf -n 3 aquilonia_fi.tab > aquilonia_fi_3inds_iter2.tab #3 random inds, 2nd iteration
# cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_ceu" {print $1}' > aquilonia_ceu.tab #6 inds
# shuf -n 3 aquilonia_fi.tab > aquilonia_all_6inds.tab; shuf -n 3 aquilonia_ceu.tab >> aquilonia_all_6inds.tab #6 inds   <<< THIS IS USED FOR THE FIRST TRY <<<
# shuf -n 3 aquilonia_fi.tab > aquilonia_all_6inds_iter2.tab; shuf -n 3 aquilonia_ceu.tab >> aquilonia_all_6inds_iter2.tab #6 random inds, 2nd iteration

# cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "Lai_1w" || $1 == "Lai_2w" || $1 == "Loa_1w" {print $1}' > aquilonia_scot.tab #3 inds
# cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CBAQ1_1w" || $1 == "CBAQ2_2w" || $1 == "CBAQ3_1w" {print $1}' > aquilonia_swi.tab #3 inds

# # F. polyctena
# cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $1}' > polyctena_ceu.tab #6 inds   <<< THIS IS USED FOR THE FIRST TRY <<<

# # F. pratensis
# cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "prat" {print $1}' > pratensis_fi.tab #7 inds   <<< THIS IS USED FOR THE FIRST TRY <<<

# # F. lugubris
# cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "lug" {print $1}' > lugubris_fi.tab #8 inds   <<< THIS IS USED FOR THE FIRST TRY <<<

### Add species information back to the files as the second column (will not re-create the files because want the same individuals for the 
### local and global Fst estimation and the group files for global Fst estimation were randomly subsampled

INPATH=/scratch/project_2001443/barriers_introgr_formica/fst_global/groupfiles
OUTPATH=/scratch/project_2001443/barriers_introgr_formica/pixy/groupfiles

# The first try
GROUPFILE=group_parentals_all_pixy.tab

paste $INPATH/rufa_all_6inds.tab <(yes "rufa" | head -n 6) > $OUTPATH/rufa_all_6inds_pixy.tab
paste $INPATH/aquilonia_all_6inds.tab <(yes "aquilonia" | head -n 6) > $OUTPATH/aquilonia_all_6inds_pixy.tab
paste $INPATH/polyctena_ceu.tab <(yes "polyctena_ceu" | head -n 6) > $OUTPATH/polyctena_ceu_6inds_pixy.tab
paste $INPATH/lugubris_fi.tab <(yes "lugubris_fi" | head -n 8) | shuf -n 6 > $OUTPATH/lugubris_fi_6inds_pixy.tab
paste $INPATH/pratensis_fi.tab <(yes "pratensis_fi" | head -n 7) | shuf -n 6 > $OUTPATH/pratensis_fi_6inds_pixy.tab

cat $OUTPATH/rufa_all_6inds_pixy.tab $OUTPATH/aquilonia_all_6inds_pixy.tab $OUTPATH/polyctena_ceu_6inds_pixy.tab \
$OUTPATH/lugubris_fi_6inds_pixy.tab $OUTPATH/pratensis_fi_6inds_pixy.tab > $OUTPATH/$GROUPFILE


# Iteration 2, if needed a replicate
# ... do with the iter2 files and new shuffling for lugubris and pratensis

