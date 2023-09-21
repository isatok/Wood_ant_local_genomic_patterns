#### Prepare group files from the sample table (same logic as with TWISST group files) ####

####NO EIKU HELVETTI, TEE NÄÄ LOPPUUN. TÄÄLLÄKI https://d197for5662m48.cloudfront.net/documents/publicationstatus/120529/preprint_pdf/cdef18babcbe3b34b8e020e6293bd445.pdf
####ON KÄYTETTY SAMANLAISESSA TILANTEESSA VCFTOOLSIA - WEIGHTED VAI MEAN SELVINNEE MYÖHEMMIN ######
#######THE MEAN AND WEIGHTED ESTIMATES VARY INSANELY LOT. WILL NOT CONTINUE THIS NOW BUT TRY TO DO THE GVCF's AND RUN FST WITH THEM########

###
### 0. Prep --------------------
###

sinteractive ...

mkdir /scratch/project_2001443/barriers_introgr_formica/fst_global
cd /scratch/project_2001443/barriers_introgr_formica/fst_global

module load biokit
bcftools query -l all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz | wc -l #101 individuals

VCFIN=/scratch/project_2001443/barriers_introgr_formica/vcf/filt/all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz #101 individuals (excluded 110, RN417, 105), singletons filtered out
FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab               #101 individuals + F.exsecta + header line

###
### 1. Compute pairwise Fst estimates --------------------
###

#within species
vcftools --gzvcf ${VCFIN} --weir-fst-pop rufa_fi.tab --weir-fst-pop rufa_ceu.tab --out rufafi_rufaceu                  #rufa_fi rufa_ceu
vcftools --gzvcf ${VCFIN} --weir-fst-pop aquilonia_fi_6inds.tab --weir-fst-pop aquilonia_ceu.tab --out aqufi_aquceu    #aquilonia_fi aquilonia_ceu (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop aquilonia_scot.tab --weir-fst-pop aquilonia_swi.tab --out aquscot_aquswi      #aquilonia_scotl aquilonia_switz

#between species
vcftools --gzvcf ${VCFIN} --weir-fst-pop rufa_all.tab --weir-fst-pop polyctena_ceu.tab --out rufa_polyctenaceu         #rufa_all polyctena_ceu #12 vs 6 - balance??!
vcftools --gzvcf ${VCFIN} --weir-fst-pop rufa_ceu.tab --weir-fst-pop polyctena_ceu.tab --out rufaceu_polyctenaceu      #rufa_ceu polyctena_ceu

vcftools --gzvcf ${VCFIN} --weir-fst-pop rufa_all.tab --weir-fst-pop aquilonia_all_6inds.tab --out rufa_polyctenaceu   #rufa_all aquilonia_all - balance??!

#polyctena_ceu aquilonia_all - balance??! ## STILL TO DO ##
#polyctena_ceu aquilonia_ceu 


###
### 2. Collect data on the mean and weighted Fst values in one file --------------------
###

            ## STILL TO DO ##

###
### 3. Info on how the sample groups were created --------------------
###

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

# All possible parental groups (by species and species further divided by geography)
# F. rufa
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa" {print $1}' > rufa_all.tab #12 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_fi" {print $1}' > rufa_fi.tab #6 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_ceu" {print $1}' > rufa_ceu.tab #6 inds

# F. aquilonia
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" {print $1}' > aquilonia_all.tab #45 inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_fi" {print $1}' > aquilonia_fi.tab #39 inds
shuf -n 6 aquilonia_fi.tab > aquilonia_fi_6inds.tab #6 random inds
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_ceu" {print $1}' > aquilonia_ceu.tab #6 inds
shuf -n 3 aquilonia_fi.tab > aquilonia_all_6inds.tab; shuf -n 3 aquilonia_ceu.tab >> aquilonia_all_6inds.tab #6 inds

cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "Lai_1w" || $1 == "Lai_2w" || $1 == "Loa_1w" {print $1}' > aquilonia_scot.tab #3 inds
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CBAQ1_1w" || $1 == "CBAQ2_2w" || $1 == "CBAQ3_1w" {print $1}' > aquilonia_swi.tab #3 inds

# F. polyctena
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $1}' > polyctena_ceu.tab #6 inds

# F. pratensis
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "prat" {print $1}' > pratensis_fi.tab #7 inds

# F. lugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "lug" {print $1}' > lugubris_fi.tab #8 inds

# F. paralugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "paral" {print $1}' > paralugubris_ceu.tab #2 inds


### END.
