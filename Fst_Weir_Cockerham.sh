#This is to calculate Fst in 20kb windows with vcftools weighted W&C estimate. If this works, I will need to do it to all species pairs. ANd/ or for the 1000 genomic blocks simulations!

cd /scratch/project_2001443/barriers_introgr_formica/fst_vcftools_windowed
sinteractive ...
module load biokit

GROUPDIR=/scratch/project_2001443/barriers_introgr_formica/fst_global/groupfiles
VCF=/scratch/project_2001443/barriers_introgr_formica/vcf/filt/all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz

#Empirical aqu-pol for samples for which we have demographic models, *20kb* non-sliding windows
vcftools --gzvcf ${VCF} --weir-fst-pop $GROUPDIR/polyctena_ceu.tab --weir-fst-pop $GROUPDIR/aquilonia_ceu.tab --fst-window-size 20000 --out polceu_aquceu_mac2_20kbW_weighted
vcftools --gzvcf ${VCF} --weir-fst-pop $GROUPDIR/polyctena_fi_4inds.tab --weir-fst-pop $GROUPDIR/aquilonia_fi_4inds.tab --fst-window-size 20000 --out polfi_aqufi_mac2_20kbW_weighted

VCFSIMFI=/scratch/project_2001443/barriers_introgr_formica/msprime/aquFI_polFI_sim_1/output_corrHeadPos.vcf.gz
VCFSIMSWI=/scratch/project_2001443/barriers_introgr_formica/msprime/aquSWI_polWSWI_sim_1/output_corrHeadPos.vcf.gz

#Simulated aqu-pol, *10kb* non-sliding windows
cd /scratch/project_2001443/barriers_introgr_formica/msprime/fst_vcftools_windowed_sim
vcftools --gzvcf ${VCFSIMFI} --weir-fst-pop group_sim_aqu.tab --weir-fst-pop group_sim_pol.tab --fst-window-size 10000 --out polfi_aqufi_10kbW_weighted_sim
vcftools --gzvcf ${VCFSIMSWI} --weir-fst-pop group_sim_aqu.tab --weir-fst-pop group_sim_pol.tab --fst-window-size 10000 --out polswi_aquswi_10kbW_weighted_sim


#Finnish polyctena samples
polyctena_fi_4inds.tab:
Att1_1w
Fis2_1w
Jar6_1w
Lok3_1w

#Finnish aquilonia samples
aquilonia_fi_4inds.tab:
CF14a_1w
CF4b_1w
CF8b_1w
Pus2_1w


#Get the data
#Empirical
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/fst_vcftools_windowed/*' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Local_Fst_vcftools'
#Simulated
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/msprime/fst_vcftools_windowed_sim/*' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/Local_Fst_vcftools'


--------------------------------------------
#Why are my Fst estimates so much lower than Pierre's & Beatrize's? Let's try without mac2 filtering and without Scaffold03.
## ---> OK this is about weighted and unweighted Fst estimates.

cd /scratch/project_2001443/barriers_introgr_formica/fst_global
sinteractive ...
mkdir filter_effect
cd filter_effect
module load biokit

VCFDIR=/scratch/project_2001443/barriers_introgr_formica/vcf/filt
GROUPDIR=/scratch/project_2001443/barriers_introgr_formica/fst_global/groupfiles

VCF_noMac=$VCFDIR/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.NoScaff00.vcf.gz
VCF_mac2=$VCFDIR/all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz
VCF_noMacNoScaff03=$VCFDIR/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.NoScaff0003.vcf.gz

vcftools --gzvcf ${VCF_noMac} --not-chr Scaffold03 --recode --stdout | bgzip > $VCF_noMacNoScaff03
bcftools index -t $VCF_noMacNoScaff03

VCF_thinned_noScaff0003_mac2=$VCFDIR/all_samples_BS.DP8.hwe.AN10.mac2.noScaff0003.thin20kb.vcf.gz

vcftools --gzvcf ${VCF_noMac} --weir-fst-pop $GROUPDIR/polyctena_ceu.tab --weir-fst-pop $GROUPDIR/aquilonia_ceu.tab --out polyctenaceu_aquiloniaceu_noMac
vcftools --gzvcf ${VCF_mac2} --weir-fst-pop $GROUPDIR/polyctena_ceu.tab --weir-fst-pop $GROUPDIR/aquilonia_ceu.tab --out polyctenaceu_aquiloniaceu_mac2   
vcftools --gzvcf ${VCF_noMacNoScaff03} --weir-fst-pop $GROUPDIR/polyctena_ceu.tab --weir-fst-pop $GROUPDIR/aquilonia_ceu.tab --out polyctenaceu_aquiloniaceu_noMacNoScaff03   
vcftools --gzvcf ${VCF_thinned_noScaff0003_mac2} --weir-fst-pop $GROUPDIR/polyctena_ceu.tab --weir-fst-pop $GROUPDIR/aquilonia_ceu.tab --out polyctenaceu_aquiloniaceu_thinned_noScaff0003_mac2   

#noMac
 #Weir and Cockerham mean Fst estimate: 0.21808
 #Weir and Cockerham weighted Fst estimate: 0.38355

#noMacNoScaff03
 #Weir and Cockerham mean Fst estimate: 0.22081
 #Weir and Cockerham weighted Fst estimate: 0.39139

#mac2
 #Weir and Cockerham mean Fst estimate: 0.22949
 #Weir and Cockerham weighted Fst estimate: 0.38822

#thinned_noScaff0003_mac2
 #Weir and Cockerham mean Fst estimate: 0.24395
 #Weir and Cockerham weighted Fst estimate: 0.4214

#### Prepare group files from the sample table (same logic as with TWISST group files) ####

###
### 0. Prep --------------------
###

sinteractive ...

#mkdir /scratch/project_2001443/barriers_introgr_formica/fst_global

module load biokit
#bcftools query -l all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz | wc -l #101 individuals

VCFIN=/scratch/project_2001443/barriers_introgr_formica/vcf/filt/all_samples.DP8.hwe.AN10.noScaff00.mac2.vcf.gz #101 individuals (excluded 110, RN417, 105), singletons filtered out
FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab               #101 individuals + F.exsecta + header line


###
### 1. Compute pairwise Fst estimates --------------------
###

cd /scratch/project_2001443/barriers_introgr_formica/fst_global

### Within each species

#rufa
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/rufa_ceu.tab --out rufafi_rufaceu #rufa_fi rufa_ceu

#polyctena
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_esw.tab --weir-fst-pop ./groupfiles/polyctena_wsw.tab --out poleastsw_polwestsw #polyctena_swi east west # <<---- new

#aquilonia
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out aqufi_aquceu    #aquilonia_fi aquilonia_ceu (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds_iter2.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out aqufi_aquceu_iter2    #aquilonia_fi aquilonia_ceu (balanced) #2nd iteration#
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_scot.tab --weir-fst-pop ./groupfiles/aquilonia_swi.tab --out aquscot_aquswi      #aquilonia_scotl aquilonia_switz
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_3inds.tab --weir-fst-pop ./groupfiles/aquilonia_scot.tab --out aqufi_aquscotl    #aquilonia_fi aquilonia_scotl (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_3inds_iter2.tab --weir-fst-pop ./groupfiles/aquilonia_scot.tab --out aqufi_aquscotl_iter2    #aquilonia_fi aquilonia_scotl (balanced) #2nd iteration#
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_3inds.tab --weir-fst-pop ./groupfiles/aquilonia_swi.tab --out aqufi_aquswitz    #aquilonia_fi aquilonia_switz (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_3inds_iter2.tab --weir-fst-pop ./groupfiles/aquilonia_swi.tab --out aqufi_aquswitz_iter2    #aquilonia_fi aquilonia_switz (balanced) #2nd iteration#
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_6_northest_fi.tab --weir-fst-pop ./groupfiles/aquilonia_6_southest_fi.tab --out aqunorthfi_aqusouthfi #aquilonia_fi north south # <<---- new

#lugubris
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_3_northest_fi.tab --weir-fst-pop ./groupfiles/lugubris_3_southest_fi.tab --out lugnorthfi_lugsouthfi #lugubris_fi north south # <<---- new

#pratensis
 #NA

### Within rufa/pol and aqu/lug clades

#rufa-pol
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/polyctena_ceu.tab --out rufa_polyctenaceu   #rufa_all polyctena_ceu (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/polyctena_ceu.tab --out rufaceu_polyctenaceu      #rufa_ceu polyctena_ceu
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/polyctena_ceu.tab --out rufafi_polyctenaceu      #rufa_fi polyctena_ceu

#aqu-lug
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_fi.tab --weir-fst-pop ./groupfiles/aquilonia_all_6inds.tab --out lugubrisfi_aquilonia   #lugubris_fi aquilonia_all (balanced) # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_fi.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out lugubrisfi_aquiloniaceu   #lugubris_fi aquilonia_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_fi.tab --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --out  lugubrisfi_aquiloniafi #lugubris_fi - aqu_fi # <<---- new


### Between rufa/pol and aqu/lug clades

#rufa-aqu
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/aquilonia_all_6inds.tab --out rufa_aquilonia   #rufa_all aquilonia_all (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds_iter2.tab --weir-fst-pop ./groupfiles/aquilonia_all_6inds_iter2.tab --out rufa_aquilonia_iter2   #rufa_all aquilonia_all (balanced) #2nd iteration#
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --out rufafi_aquiloniafi #rufa_fi - aqu_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out rufaceu_aquiloniaceu #rufa_ceu - aqu_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out rufafi_aquiloniaceu #rufa_fi - aqu_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --weir-fst-pop ./groupfiles/rufa_ceu.tab --out aquiloniafi_rufaceu #aqu_fi - rufa_ceu # <<---- new

#pol-aqu
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_all_6inds.tab --out polyctenaceu_aquilonia   #polyctena_ceu aquilonia_all (balanced)
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out polyctenaceu_aquiloniaceu   #polyctena_ceu aquilonia_ceu 
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --out  polyctenaceu_aquiloniafi #pol_ceu - aqu_fi # <<---- new

#rufa-lug
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out  rufafi_lugubrisfi #rufa_fi lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out  rufaceu_lugubrisfi #rufa_ceu lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out rufa_lugubrisfi #rufa_all lug_fi # <<---- new
 
#pol-lug
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out polyctenaceu_lugubrisfi #pol_ceu - lug_fi # <<---- new
 
### All species compared to pratensis

#aqu-prat
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out aquiloniafi_pratensisfi #aqu_fi - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out aquiloniaceu_pratensisfi #aqu_ceu - prat_fi # <<---- new
 
#lug-prat
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_fi.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out lugubrisfi_pratensisfi #lug_fi prat_fi # <<---- new
 
#pol-prat
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out polyctenaceu_pratensisfi #pol_ceu prat_fi # <<---- new
 
#rufa-prat
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufaceu_pratensisfi #rufa_ceu - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufafi_pratensisfi #rufa_fi - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufa_pratensisfi #rufa_all - prat_fi # <<---- new

#run these new ones
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_esw.tab --weir-fst-pop ./groupfiles/polyctena_wsw.tab --out poleastsw_polwestsw #polyctena_swi east west # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_6_northest_fi.tab --weir-fst-pop ./groupfiles/aquilonia_6_southest_fi.tab --out aqunorthfi_aqusouthfi #aquilonia_fi north south # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_3_northest_fi.tab --weir-fst-pop ./groupfiles/lugubris_3_southest_fi.tab --out lugnorthfi_lugsouthfi #lugubris_fi north south # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --out rufafi_aquiloniafi #rufa_fi - aqu_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out rufaceu_aquiloniaceu #rufa_ceu - aqu_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --out rufafi_aquiloniaceu #rufa_fi - aqu_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --weir-fst-pop ./groupfiles/rufa_ceu.tab --out aquiloniafi_rufaceu #aqu_fi - rufa_ceu # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --out  polyctenaceu_aquiloniafi #pol_ceu - aqu_fi # <<---- new
poleastsw_polwestsw
aqunorthfi_aqusouthfi


# check the results this far first - is everything ok?
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out  rufafi_lugubrisfi #rufa_fi lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out  rufaceu_lugubrisfi #rufa_ceu lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out rufa_lugubrisfi #rufa_all lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/lugubris_fi.tab --out polyctenaceu_lugubrisfi #pol_ceu - lug_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_fi_6inds.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out aquiloniafi_pratensisfi #aqu_fi - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/aquilonia_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out aquiloniaceu_pratensisfi #aqu_ceu - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/lugubris_fi.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out lugubrisfi_pratensisfi #lug_fi prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/polyctena_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out polyctenaceu_pratensisfi #pol_ceu prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_ceu.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufaceu_pratensisfi #rufa_ceu - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_fi.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufafi_pratensisfi #rufa_fi - prat_fi # <<---- new
vcftools --gzvcf ${VCFIN} --weir-fst-pop ./groupfiles/rufa_all_6inds.tab --weir-fst-pop ./groupfiles/pratensis_fi.tab --out rufa_pratensisfi #rufa_all - prat_fi # <<---- new


###
### 2. Collect data on the mean and weighted Fst values in one file --------------------
###

cd /scratch/project_2001443/barriers_introgr_formica/fst_global/

output_file="fst_results.txt" 
column_names=("Weir and Cockerham mean Fst" "Weir and Cockerham weighted Fst")

> "$output_file"                                                 #emtpy the output file if it exists
echo -e "Filename\t${column_names[@]}" >> "$output_file"         #add a header to the output file

for file in *.log; do             #iterate through all *.log files
filename="${file%.log}"           #extract the filename excluding ".log"
values=()                         #initialise an array for the fst values
while IFS= read -r line; do       #read each line in the file
 for col_name in "${column_names[@]}"; do
        if [[ $line == "$col_name"* ]]; then #these check that the line starts with the defined phrases (same as column_names)
 value=$(echo "$line" | grep -oP '[0-9]+\.[0-9]+')
 values+=("$value")
        fi
     done
done < "$file"

echo -e "$filename\t${values[@]}" >> "$output_file" #print results to the output file
done

scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/fst_global/fst_results.txt' \
'/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/fst_global'

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

cut -f1,4,5 $FULLSAMPLE | sort -k3 | awk '$2 == "aqu_fi" {print $0}' | head -n 6 | awk '{print $1}' > aquilonia_6_southest_fi.tab #6 most southern Finnish aqu (>6degr lat diff to northest)
cut -f1,4,5 $FULLSAMPLE | sort -k3 | awk '$2 == "aqu_fi" {print $0}' | tail -n 6 | awk '{print $1}' > aquilonia_6_northest_fi.tab #6 most northern Finnish aqu

# F. polyctena
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "pol" {print $1}' > polyctena_ceu.tab #6 inds

cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CBCH1_1w" || $1 == "CBCH2_2w" || $1 == "CBCH3_1w" {print $1}' > polyctena_esw.tab #3 inds, eastern switz (ca 3degr lon diff to wsw)
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' | awk '$1 == "CAGa_1w" || $1 == "NAZa_1w" || $1 == "VDa_1w" {print $1}' > polyctena_wsw.tab #3 inds, western switz

# F. pratensis
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "prat" {print $1}' > pratensis_fi.tab #7 inds

# F. lugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "lug" {print $1}' > lugubris_fi.tab #8 inds

cut -f1,4,5 $FULLSAMPLE | sort -k3 | awk '$2 == "lug_fi" {print $0}' | head -n 3 | awk '{print $1}' > lugubris_3_southest_fi.tab #3 most southern Finnish lug (min >2degr lat diff to northest)
cut -f1,4,5 $FULLSAMPLE | sort -k3 | awk '$2 == "lug_fi" {print $0}' | tail -n 3 | awk '{print $1}' > lugubris_3_northest_fi.tab #3 most northern Finnish lug

# F. paralugubris
cut -f1,3 $FULLSAMPLE | grep -v 'vcf_id' | awk '$2 == "paral" {print $1}' > paralugubris_ceu.tab #2 inds


### END.
