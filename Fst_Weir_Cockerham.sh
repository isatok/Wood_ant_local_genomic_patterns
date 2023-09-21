#### Prepare group files from the sample table (same logic as with TWISST group files) ####

mkdir /scratch/project_2001443/barriers_introgr_formica/fst_global

cd /scratch/project_2001443/barriers_introgr_formica/fst_global

#Use singleton-filtered full vcf where samples 110, RN417 and 105 are filtered out
module load biokit
bcftools query -l all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz | wc -l #101 individuals

VCFIN=/scratch/project_2001443/barriers_introgr_formica/vcf/filt/all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz #101 individuals
FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab               #101 individuals + F.exsecta + header line

#Fst pairs
vcftools --gzvcf ${VCFIN} --weir-fst-pop rufa_all.tab --weir-fst-pop polyctena_ceu.tab --out rufa_polyctena #rufa-polyctena all
#rufa_ceu polyctena_ceu
#...


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
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa" {print $1}' > rufa_all.tab
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_fi" {print $1}' > rufa_fi.tab
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "rufa_ceu" {print $1}' > rufa_ceu.tab

# F. aquilonia
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu" {print $1}' > aquilonia_all.tab
cut -f1,4 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "aqu_fi" {print $1}' > aquilonia_fi.tab
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "Lai_1w" || $1 == "Lai_2w" || $1 == "Loa_1w" {print $1}' > aquilonia_scot.tab
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CBAQ1_1w" || $1 == "CBAQ2_2w" || $1 == "CBAQ3_1w" {print $1}' > aquilonia_swi.tab

# F. polyctena
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $1}' > polyctena_ceu.tab

# F. pratensis
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "prat" {print $1}' > pratensis_fi.tab

# F. lugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "lug" {print $1}' > lugubris_fi.tab

# F. paralugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "paral" {print $1}' > paralugubris_ceu.tab
