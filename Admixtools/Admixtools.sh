##env.yml
#name: admixtools_env
#channels:
#  - conda-forge
#  - bioconda
#  - defaults
#dependencies:
#  - admixtools

#Use Admixtools:
  export PATH="/projappl/project_2001443/admixtools_env/bin:$PATH" 

# Is this important?
#  "([ INFO ] Creating wrappers 
#FATAL:   exec /usr/bin/test failed: input/output error)"



#Filter out Scaffold00
sinteractive --account project_2001443 --mem 4000

VCFIN=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.vcf.gz
VCFOUT=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.NoScaff00.vcf.gz
VCFPATH=/scratch/project_2001443/barriers_introgr_formica/vcf/filt

cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt/
module load biokit

vcftools --gzvcf $VCFIN --not-chr Scaffold00 --remove-indv 105-FaquH --remove-indv 110-FaquH --recode --stdout | bgzip > $VCFOUT #Filtering out Scaff00 is unnecessary; it's not been SNP-called
bcftools index -t $VCFOUT
bcftools index -n $VCFOUT #708783 SNPs.

#So now we have a VCF file for Admixtools that has all SNPs in all "proper" chromosomes. No thinning, no MAC filtering. Samples 105-FaquH and 110-FaquH filtered out.

#Convert the VCF to eigenstrat format using Joana Meier's script

cd /scratch/project_2001443/barriers_introgr_formica/admixtools
mv $VCFPATH/$VCFOUT .

sh convertVCFtoEigenstrat.sh all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.NoScaff00.vcf.gz

#I set the recombination rate to 5 (cM/Mb), since the average rho from Nouhaud 2022 is 0.0049 (I assume cM/bp??),
#which is in line with known eusocial Hymenopteran recombination rates (5-7; https://doi.org/10.1073/pnas.1208094109),
#and changed renameScaff="TRUE", and added Admixtools PATH to the shell script.
# #of SNPs left: numsnps output: 695502

#Then, in R:
#Create a directory for this analysis (in terminal):
mkdir /users/satokan1/privatemodules/admixr

#edit the .Renviron file as guided in https://bodkan.net/admixr/articles/tutorial.html (in console):
usethis::edit_r_environ()
PATH="/projappl/project_2001443/admixtools_env/bin:$PATH" 

#download admixr (in terminal):
cd /projappl/project_2001443
mkdir project_rpackages_4.2.1 

.libPaths(c("/projappl/project_2001443/project_rpackages_4.2.1", .libPaths()))
libpath <- .libPaths()[1]
#install.packages("admixr", lib = libpath) #does not work - says the admixr is not compatible with my R version, but the reason may be more complicated (googling)
library(devtools)


library(admixr)
library(tidyverse)

# set data prefix
data_prefix <- "./admixtools/sparrows"
# read in data
snps <- eigenstrat(data_prefix)
# count SNPs
snp_info <- count_snps(snps)


#... look more from https://speciationgenomics.github.io/ADMIXTOOLS_admixr/. Need to modify the pop names at least! And download & make the program work in Puhti R.

