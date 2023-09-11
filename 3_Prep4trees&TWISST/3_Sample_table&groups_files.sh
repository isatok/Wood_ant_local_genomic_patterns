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

###Read in the sample list and create a data frame
samples <- read.csv("sample.list.csv", header=F)
samples2 <- as.data.frame(t(samples[1,]))

###Create vectors of samples which belong to each species or group

###Lists by species
prat <- c("117-Fprat", "120-Fprat", "122-Fprat", "123-Fprat", "2-Fprat", "3-Fprat", "6-Fprat")
lug <- c("1-Fprat", "108-Flug", "11-Flug", "113-Flug", "115-Flug", "124-Flug", "56-Flug", "88-Flug")
rufa <- c("102-Frufa", "12-Frufa", "16-Frufa", "26-Frufa", "37-Frufa", "72-Frufa")
pol <- c("CAGa_1w", "CBCH1_1w", "CBCH2_2w", "CBCH3_1w", "NAZa_1w", "VDa_1w")
aqu <- c("10-Fpol", "103-Faqu", "104-Faqu", "106-Faqu", "107-Faqu", "109-Faqu", "114-FaquH", "13-Faqu", "17-Faqu", "22-Faqu", "24-Faqu", "25-Faqu", "27-Fpol", "30-Faqu", "31-Faqu", "32-Frufa", "38-Fpol", "42-Faqu", "44-FaquH", "47-FaquH", "48-Fpol", "50-FaquH", "51-FaquH", "58-Fpol", "62-FaquH", "66-Fpol", "67-Faqu", "70-Faqu", "76-FaquH", "77-Fpol", "83-Faqu", "9-Faqu", "91-FpolH", "93-Faqu", "99-Faqu", "CF14a_1w", "CF4b_1w", "CF8b_1w", "Pus2_1w", "CBAQ1_1w", "CBAQ2_2w", "CBAQ3_1w", "Lai_1w", #"Lai_2w", "Loa_1w")
apr_hybr <- c("105-FaquH", "14-Fpol", "36-FpolH", "43-Fpol", "52-Fpol", "53-Fpol", "54-Frufa", "55-FaquH", "59-FpolH", "65-Frufa", "Att1_1w", "Fis2_1w", "Jar6_1w", "Lok3_1w")
l_hybr <- c("119-Flug", "125-Flug", "126-Flug")

###Lists by species and by geography
pol_sw <- c("CAGa_1w", "CBCH1_1w", "CBCH2_2w", "CBCH3_1w", "NAZa_1w", "VDa_1w")
aqu_swsc <- c("CBAQ1_1w", "CBAQ2_2w", "CBAQ3_1w", "Lai_1w", "Lai_2w", "Loa_1w")
aqu_fi <- c("10-Fpol", "103-Faqu", "104-Faqu", "106-Faqu", "107-Faqu", "109-Faqu", "114-FaquH", "13-Faqu", "17-Faqu", "22-Faqu", "24-Faqu", "25-Faqu", "27-Fpol", "30-Faqu", "31-Faqu", "32-Frufa", "38-Fpol", "42-Faqu", "44-FaquH", "47-FaquH", "48-Fpol", "50-FaquH", "51-FaquH", "58-Fpol", "62-FaquH", "66-Fpol", "67-Faqu", "70-Faqu", "76-FaquH", "77-Fpol", "83-Faqu", "9-Faqu", "91-FpolH", "93-Faqu", "99-Faqu", "CF14a_1w", "CF4b_1w", "CF8b_1w", "Pus2_1w")
prat_fi <- c("117-Fprat", "120-Fprat", "122-Fprat", "123-Fprat", "2-Fprat", "3-Fprat", "6-Fprat")
lug_fi <- c("1-Fprat", "108-Flug", "11-Flug", "113-Flug", "115-Flug", "124-Flug", "56-Flug", "88-Flug")
rufa_fi <- c("102-Frufa", "12-Frufa", "16-Frufa", "26-Frufa", "37-Frufa", "72-Frufa")
apr_hybr_fi <- c("105-FaquH", "14-Fpol", "36-FpolH", "43-Fpol", "52-Fpol", "53-Fpol", "54-Frufa", "55-FaquH", "59-FpolH", "65-Frufa", "Att1_1w", "Fis2_1w", "Jar6_1w", "Lok3_1w")
l_hybr_fi <- c("119-Flug", "125-Flug", "126-Flug")

###Build the sample table
colnames(samples2) <- c("vcfID")
samples2$species[samples2$vcfID %in% prat] <- "pratensis"
samples2$species[samples2$vcfID %in% lug] <- "lugubris"
samples2$species[samples2$vcfID %in% rufa] <- "rufa"
samples2$species[samples2$vcfID %in% pol] <- "polyctena"
samples2$species[samples2$vcfID %in% aqu] <- "aquilonia"
samples2$species[samples2$vcfID %in% apr_hybr] <- "aqupolrufa_hybrids"
samples2$species[samples2$vcfID %in% l_hybr] <- "lugubris_hybrids"

samples2$species_geog[samples2$vcfID %in% pol_sw] <- "polyctena_sw"
samples2$species_geog[samples2$vcfID %in% aqu_swsc] <- "aquilonia_swsc"
samples2$species_geog[samples2$vcfID %in% prat_fi] <- "pratensis_fi"
samples2$species_geog[samples2$vcfID %in% lug_fi] <- "lugubris_fi"
samples2$species_geog[samples2$vcfID %in% rufa_fi] <- "rufa_fi"
samples2$species_geog[samples2$vcfID %in% aqu_fi] <- "aquilonia_fi"
samples2$species_geog[samples2$vcfID %in% apr_hybr_fi] <- "aqupolrufa_hybrids_fi"
samples2$species_geog[samples2$vcfID %in% l_hybr_fi] <- "lugubris_hybrids_fi"

###Export into "sample_table.tab" to be used in making TWISST groups files
#write.table(samples2, file="sample_table.tab", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


####
#### 2. Create groups files for TWISST -------------------------------------------------------------------------------------------
####

#### Prepare group files from the sample table ####

###possible groups ($2 & $3)
#aquilonia	aquilonia_fi
#aquilonia	aquilonia_swsc
#polyctena	polyctena_sw
#rufa    rufa_fi
#aqupolrufa_hybrids	aqupolrufa_hybrids_fi
#lugubris        lugubris_fi
#lugubris_hybrids        lugubris_hybrids_fi
#pratensis	pratensis_fi
#exsecta  exsecta

cd $SCRATCH/analysis/twisst/
FULLSAMPLE=/scratch/project_2001443/vcf/ind_lists/sample_table.tab

##### A) AFTER CREATING ANY ONE OF THE GROUPTMP, GO TO B) TO ADD A&B TO TIPS #####

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

cat $GROUPTMP > tmp
perl -npe 's/\t/_A\t/' tmp > tmpA
perl -npe 's/\t/_B\t/' tmp > tmpB
cat tmpA tmpB | sort > $GROUPFILE
rm tmp tmpA tmpB $GROUPTMP #remove exsecta B & letter A when present

