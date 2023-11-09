### 1. Install AdmixTools

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
      # or use admixtools2_env if the above doesn't work; the first one had the error message below
      #  "([ INFO ] Creating wrappers 
      #FATAL:   exec /usr/bin/test failed: input/output error)"


### 2. Filter VCF: filter out Scaffold00 & collaborative/ poor quality samples
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


### 3. Convert the VCF to eigenstrat format using Joana Meier's script

cd /scratch/project_2001443/barriers_introgr_formica/admixtools
mv $VCFPATH/$VCFOUT .

sh convertVCFtoEigenstrat.sh all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.NoScaff00.vcf.gz

#I set the recombination rate to 5 (cM/Mb), since the average rho from Nouhaud 2022 is 0.0049 (I assume cM/bp??),
#which is in line with known eusocial Hymenopteran recombination rates (5-7; https://doi.org/10.1073/pnas.1208094109),
#and changed renameScaff="TRUE", and added Admixtools PATH to the shell script.
# #of SNPs left: numsnps output: 695502

#Then, in R:
#Create a directory for this analysis (in terminal):  ###IS this step needed?
mkdir /users/satokan1/privatemodules/admixr ###?


### 4. Install Admixr and make it work together with AdmixTools

# a) Edit the .Renviron file as guided in https://bodkan.net/admixr/articles/tutorial.html (in console),
#    to provide R information on where to find AdmixTools:
usethis::edit_r_environ()
PATH="/projappl/project_2001443/admixtools_env/bin" # Does not work, or is sufficient by itself?

# b) Needed to fix the PATH issue by giving (also?) a "PATH" variable in R (console):
PATH="/projappl/project_2001443/admixtools_env/bin:$PATH" 

# c) Download admixr (in terminal):
cd /projappl/project_2001443
mkdir project_rpackages_4.3.0 #need this R version to make it work; 4.2 does not work.

# d) In R, add the path to Admixr package to libPath:
.libPaths(c("/projappl/project_2001443/project_rpackages_4.3.0", .libPaths()))
libpath <- .libPaths()[1]

# e) install Admixr
install.packages("admixr", lib = libpath) 

### 5. Run Admixr! 

# Load libraries. Need to give at least the PATH variable again when starting a new session.
library(admixr)
library(tidyverse)

# set data prefix
data_prefix <- "./admixtools/sparrows"
# read in data
snps <- eigenstrat(data_prefix)
# count SNPs
snp_info <- count_snps(snps)


#... look more from https://speciationgenomics.github.io/ADMIXTOOLS_admixr/. Need to modify the pop names at least!

