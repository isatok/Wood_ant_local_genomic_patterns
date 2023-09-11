#parceVCF.py script:
/scratch/project_2001443/analysis/genomics_simon/genomics_general/VCF_processing/parseVCF.py

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/

#Phased vcf, with F. exsecta added in:
FULLVCF=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/phased_with_outgroup_gtfix.vcf.gz

#Sample table (created with R script below)
sample_table.tab

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

#create a ploidy table
cut -f1 sample_table.tab | grep -v 'vcf_id' | sed 's/$/ 2/' > ploidy.tab
nano ploidy.tab #Change manually F. exsecta to "1"
PLOIDYTAB=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/ploidy.tab

#make a geno file
PARSEPY=/scratch/project_2001443/analysis/genomics_simon/genomics_general/VCF_processing/parseVCF.py

python $PARSEPY -i $FULLVCF --skipIndels --excludeDuplicates --ploidyFile $PLOIDYTAB | bgzip > phased.exsecta.geno.gz
