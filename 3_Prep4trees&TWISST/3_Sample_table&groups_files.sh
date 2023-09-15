####
#### 1. Create a sample table -------------------------------------------------------------------------------------------
####

#### Create a sample list from phased VCF ####

#parceVCF.py script:
/scratch/project_2001443/analysis/genomics_simon/genomics_general/VCF_processing/parseVCF.py

#Phased vcf:

PHASEDVCF=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/phased_with_outgroup_gtfix.vcf.gz

#create a sample list from the phased vcf file:
cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
gunzip -c $PHASEDVCF | grep 'CHROM' | perl -npe 's/\#CHROM.+FORMAT//' | \
sed 's/      /,/2g' -> sample.list.csv #to make sed recongize "tab", insert it in console with ctrl+v+[Tab]

#move it to local computer to add species and geographical information
scp satokan1@puhti.csc.fi:'/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample.list.csv' '/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/twisst_autumn2023'


#### in R - build a species-informative table from the list ####

  library(tidyverse)
  setwd("/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/twisst_autumn2023")
# Read in the sample list created from the analysed vcf (Phased vcf, with F. exsecta), and create a data frame
  samples <- read.csv("sample.list.csv", header=F, strip.white = T)
  samples2 <- as.data.frame(t(samples[1,]))
# Read in the sample metadata. (eg species) from the sample table used for all analyses
  analyslist <- read.csv("/Users/inaukkar/Library/CloudStorage/OneDrive-UniversityofHelsinki/PhD/4_formica_local_genomics/analysis_sample_list.txt", header=T, sep="\t")
# Join the sample info tables and rename vcf_id
  samples3 <- left_join(samples2, analyslist, by = c("1"="vcf_id"))
  colnames(samples3)[1] <- "vcf_id"
# Fill in information for exsecta
  samples3$id[samples3$vcf_id %in% "Fexs"] <- "Fexs"
  samples3$group[samples3$vcf_id %in% "Fexs"] <- "exsecta"
  samples3$group_geo[samples3$vcf_id %in% "Fexs"] <- "exsecta"
# Export into "sample_table.tab" to be used in making TWISST groups files
  write.table(samples3, file="sample_table.tab", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


####
#### 2. Create groups files for TWISST -------------------------------------------------------------------------------------------
####

#### Prepare group files from the sample table ####

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
FULLSAMPLE=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/sample_table.tab

cut -f3 $FULLSAMPLE | sort | uniq   #check which groups are possible (w-out geogr. information)
# aqu
# aquxpolrufa
# exsecta
# lug
# lugxunk
# paral
# pol
# prat
# rufa

cut -f4 $FULLSAMPLE | sort | uniq   #check which groups are possible (WITH geogr. information)
# adm1
# adm2
# aqu_ceu
# aqu_fi
# exsecta
# lug_ceu
# lug_fi
# paral_ceu
# pol_ceu
# prat_fi
# rufa_ceu
# rufa_fi



##### A) AFTER CREATING ANY ONE OF THE GROUPTMP, GO TO B) TO ADD A&B TO TIPS #####

### all samples, to be used as the popsfile in freq.py
GROUPFILE=all_samples_pops.txt
cut -f1 $FULLSAMPLE | grep -v 'vcf_id' > /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###rufa, aquilonia, polyctena, admixed1 -> see the general pattern? think about this...


##### THESE GROUPFILES BELOW ARE FROM THE OLDER 2022 ANALYSIS RUN. CAN USE IF NEEDED. 

###non-admixed
GROUPFILE=group_parentals.tab
cut -f1,2 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia" || \
$2 == "polyctena" || $2 == "rufa" || $2 == "lugubris" || $2 == "pratensis" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###non-admixed Finnish
GROUPFILE=group_parentals_fi.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "pratensis_fi" || \
$2 == "aquilonia_fi" || $2 == "rufa_fi" || $2 == "lugubris_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp


###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_1"
GROUPFILE=group_hybrid_1.tab
cut -f1,2 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia" || \
$2 == "polyctena" || $2 == "rufa" || $2 == "aqupolrufa_hybrids" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_2"
GROUPFILE=group_hybrid_2.tab
cut -f1,2 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia" || \
$2 == "polyctena" || $2 == "aqupolrufa_hybrids" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_3"
GROUPFILE=group_hybrid_3.tab
cut -f1,2 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia" || \
$2 == "rufa" || $2 == "aqupolrufa_hybrids" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_4"
GROUPFILE=group_hybrid_4.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia_swsc" || \
$2 == "aquilonia_fi" || $2 == "rufa_fi" || $2 == "aqupolrufa_hybrids_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_5"
GROUPFILE=group_hybrid_5.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia_swsc" || \
$2 == "aquilonia_fi" || $2 == "rufa_fi" || $2 == "aqupolrufa_hybrids_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_6" DONE
GROUPFILE=group_hybrid_6.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia_swsc" || \
$2 == "polyctena_sw" || $2 == "aqupolrufa_hybrids_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_7"
GROUPFILE=group_hybrid_7.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia_fi" || \
$2 == "rufa_fi" || $2 == "aqupolrufa_hybrids_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_8"

###reduced parentals, exsecta, aqupolrufa_hybrids: "hybrid_9"
GROUPFILE=group_hybrid_9.tab
cut -f1,3 $FULLSAMPLE | grep -v 'vcfID' | awk '$2 == "aquilonia_swsc" || \
$2 == "aquilonia_fi" || $2 == "polyctena_sw" || $2 == "rufa_fi" || $2 == "aqupolrufa_hybrids_fi" || $2 == "exsecta" {print $0}' > $SCRATCH/analysis/twisst/$GROUPFILE.tmp
GROUPTMP=$GROUPFILE.tmp


##### B) - GENERAL #####

### all samples, to be used as the popsfile in freq.py #ACTUALLY THIS IS NOT NEEDED BUT NOW I CREATE THE PLOIDYAB THROUGH THIS#

cat $GROUPTMP > tmp
perl -npe 's/$/_A/' tmp > tmpA
perl -npe 's/$/_B/' tmp > tmpB
cat tmpA tmpB | sort > tmpC
cat tmpA tmpB | sort > tmpD
paste tmpC tmpD > $GROUPFILE
#THEN manually move the outgroup (Fexs) as the last one on the list, as the last one is interpreted as the outgroup.

### also make a ploidyAB.tab to be used by freq.py
cut -f1 all_samples_pops.txt | perl -npe 's/$/ 1/' > ploidyAB.tab 


### Groups files for TWISST

cat $GROUPTMP > tmp
perl -npe 's/\t/_A\t/' tmp > tmpA
perl -npe 's/\t/_B\t/' tmp > tmpB
cat tmpA tmpB | sort > $GROUPFILE
rm tmp tmpA tmpB $GROUPTMP #remove exsecta B & letter A when present

