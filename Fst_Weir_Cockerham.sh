#### Prepare group files from the sample table (same logic as with TWISST group files) ####

mkdir scratch/project_2001443/barriers_introgr_formica/fst
mkdir scratch/project_2001443/barriers_introgr_formica/fst/global
mkdir scratch/project_2001443/barriers_introgr_formica/fst/local



cd scratch/project_2001443/barriers_introgr_formica/fst/global

#Use singleton-filtered full vcf where samples 110, 105, and 54 are filtered out (or not used)
VCFIN= XXXXXXXX
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

# All possible parental groups (by geography or by species only)
# F. rufa
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa" {print $0}' > rufa_all.tab
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_fi" {print $0}' > rufa_fi.tab
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_ceu" {print $0}' > rufa_ceu.tab

# F. aquilonia
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" {print $0}' > aquilonia_all.tab
aqu fi
aqu ceu
aqu swi
aqu scotl

# F. polyctena
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $0}' > polyctena_ceu.tab

# F. pratensis
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "prat" {print $0}' > pratensis_fi.tab

# F. lugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "lug" {print $0}' > lugubris_fi.tab

# F. paralugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "paral" {print $0}' > paralugubris_ceu.tab
