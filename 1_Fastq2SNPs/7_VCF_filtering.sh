# 6_VCF_filtering.sh


## Since different parts of the pipeline require different resources, the master script here is split in three scripts

######################################################## Beginning script 1

#!/bin/bash -l
#SBATCH -J filt1
#SBATCH --account=project_2001443
#SBATCH -t 48:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=12G

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"
module load biokit

cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt
TMPDIR=/scratch/project_2001443/barriers_introgr_formica/tmp

##
## 0. Pre-filtering: normalisation, indels, non-SNPs, read imbalance, decomposition ---------------
##

vt normalize -n -r ../../../reference_genome/Formica_hybrid_v1_wFhyb_Sapis.fa ../raw/all_samples.vcf.gz | bgzip -c > all_samples.normalized.vcf.gz

bcftools filter --threads 8 -Oz -s+ --SnpGap 2 all_samples.normalized.vcf.gz > all_samples.normalized.SnpGap_2.vcf.gz && \

bcftools filter --threads 8 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ all_samples.normalized.SnpGap_2.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 8 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 & SAF>=1 && SAR>=1' all_samples.normalized.SnpGap_2.NonSNP.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 8 -O z -f PASS all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 8 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir $TMPDIR --max-mem 4G -O z > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \

bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz




##
## 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
##

bcftools filter --threads 8 --include 'QUAL >= 30 && TYPE="snp"' -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 8 --min-alleles 2 --max-alleles 2 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz




##
## 2. Prep for max DP filtering -------------------------------------------------------------------
##

### a. Correct field IDs (MNV issue)

# Extract header from VCF
bcftools view -h all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

# Check fraction missing data
vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz --missing-indv --out all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.MISSING_DATA

echo "Average missing data over all individuals:"
cut -f5 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.MISSING_DATA.imiss | grep -v 'F_' | awk '{sum+=$1;} END{print sum/NR;}'


### b. Compute mean DP per individual
vcftools --depth --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz --out all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader


### c. Compute individual DP thresholds
module load r-env
Rscript minMaxDP.R

# minMaxDP.R script copied below, for the record

      # dat = read.table("all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth", h = T)
      # dat$MIN_DEPTH = round((dat$MEAN_DEPTH-.5)/2) # min depth is half of mean depth
      # dat$MAX_DEPTH = round((dat$MEAN_DEPTH)*2) # max depth is twice mean depth
      # write.table(dat, "all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth.MinMax", col.names = T, row.names = F, quote = F, sep = "\t")


### d. Create results folder
mkdir logs
mkdir individualDP


######################################################## End script 1




##
## 3. Filtering on sequencing depth ---------------------------------------------------------------
##

#create a list of individuals remaining in the vcf
cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt
vcf-query -l all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > ind.list

######################################################## Beginning script 2

#!/bin/bash -l
#SBATCH -J splitDP
#SBATCH --account=project_2001443
#SBATCH -o logs/splitDP_%a.out
#SBATCH -e logs/splitDP_%a.err
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-103
#SBATCH --ntasks 1
#SBATCH --mem=2G

module load biokit
cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt

# Get file name and sample ID
ind=$(sed -n "$SLURM_ARRAY_TASK_ID"p ind.list)
mkdir individualDP/${ind}
cd individualDP/${ind}

# Extract single individual from VCF and apply DP filters

dpmin=2 # hard-coded, we will apply the same min filter accross all individuals afterwards
dpmax=$(grep -w "$ind" ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth.MinMax | cut -f5)

echo "Processing individual $ind with $dpmin < DP <= $dpmax"
bcftools view -s ${ind} -Ou ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz | bcftools filter -e "FORMAT/DP < ${dpmin} | FORMAT/DP >= ${dpmax}" --set-GTs . -Oz > ${ind}.vcf.gz &&
bcftools index -t ${ind}.vcf.gz

######################################################## End script 2




##
## 4. HWE testing ---------------------------------------------------------------------------------
##

######################################################## Beginning script 3

#!/bin/bash -l
#SBATCH -J filt3
#SBATCH --account=project_2001443
#SBATCH -t 48:00:00
#SBATCH -p small
#SBATCH --ntasks 4
#SBATCH --mem=8G


cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt/individualDP
module load biokit
module load r-env


### Combine all samples
ls */*gz > sample.list
bcftools merge -l sample.list -Oz --threads 4 > ../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz


### Test HWE in the samples
cd ..
vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz --hardy --out SAMPLES_hwe

### Extract sites with heterozygote excess, in R
Rscript hwe_filtering.R

# hwe_filtering.R copied below, for the record

      # # get the filename to read the hwe test
      # filename = "SAMPLES_hwe" # file name of reference table
      # # print the command line options
      # print(paste("file name with hwe output:", filename))
      #
      # # Read the HW results
      # hwvalues <- read.table(paste(filename,".hwe",sep=""), header=T, stringsAsFactors = F)
      # # Get distribution of the p-values
      # print(paste("he excess", sum(hwvalues$P_HET_EXCESS<0.01)))
      # print(paste("he deficit", sum(hwvalues$P_HET_DEFICIT<0.01)))
      # print(paste("he overall", sum(hwvalues$P_HWE<0.01)))
      #
      # # remove only the sites with excess of heterozygotes
      # indextoremove = which(hwvalues$P_HET_EXCESS<0.01)
      # # Create BED file with the sites that fail the HW test for excess of heterozygotes
      # # BED file has three entries:
      # # chromosome, position-1, position (it is position-1 because it is assumed that the first base is 0)
      # position = hwvalues[indextoremove, 2]
      # bedmatrix = matrix(c(hwvalues[indextoremove, 1], format(position-1, scientific = F), format(position, scientific = F)), ncol=3, byrow=F)
      # write("#Sites with heterozygosity excess", file=paste("nohwe_excess_",filename,".bed",sep=""))
      # write(t(bedmatrix), file=paste("nohwe_excess_",filename,".bed",sep=""), ncolumns = 3, append=T)


### Remove these sites in the VCF

vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz --exclude-bed nohwe_excess_SAMPLES_hwe.bed --recode --recode-INFO-all --stdout | bgzip >  all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz




##
## 5. Min DP filtering  ---------------------------------------------------------------------------
##

# Individual genotypes not reaching the minimal cutoff dpmin are tagged as missing (--set-GTs .)

# DP 5
dpmin=5
bcftools filter --threads 4 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz

# DP 6
dpmin=6
bcftools filter --threads 4 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz

# DP 7
dpmin=7
bcftools filter --threads 4 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz

# DP 8
dpmin=8
bcftools filter --threads 4 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz


##
## 6. Filter based on missing data ----------------------------------------------------------------
##


### FILTERING BASED ON ALLELIC NUMBER (AN)

# 103 diploid samples = 206 alleles expected 
# 10% missing data = 185,4 = 186 alleles min.

# DP5
inputfile=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP5.hwe.vcf.gz
bcftools view --threads 4 -Oz -i 'AN >= 186' ${inputfile} > ${inputfile%.vcf.gz}.AN10percMiss.vcf.gz

# DP6
inputfile=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP6.hwe.vcf.gz
bcftools view --threads 4 -Oz -i 'AN >= 186' ${inputfile} > ${inputfile%.vcf.gz}.AN10percMiss.vcf.gz

# DP7
inputfile=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP7.hwe.vcf.gz
bcftools view --threads 4 -Oz -i 'AN >= 186' ${inputfile} > ${inputfile%.vcf.gz}.AN10percMiss.vcf.gz

# DP8
inputfile=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.vcf.gz
bcftools view --threads 4 -Oz -i 'AN >= 186' ${inputfile} > ${inputfile%.vcf.gz}.AN10percMiss.vcf.gz


# Index VCFs, count number of sites left overall & percentage of missing data
for i in *percMiss.vcf.gz ; do
  echo $i ;
  bcftools index -ft $i
  bcftools index -n $i
  vcftools --gzvcf ${i} --missing-indv --out ${i%.vcf.gz}
done

######################################################## End script 3



############################## Filter again AN after dropping out individuals ### Begin script 4 ########################

#After filtering for depth (DP), we excluded individuals with more than 50% missing data before we continued to filtering by allelic number. 

#We identified seven individuals with more than 50% missing data (RN418:74%, 105-FaquH:81%, 54-Frufa:63%, RN421:63%, RN425:61%, RN426:56%, RN422:56%)
#We excluded these alongside with three individual F. rufa group genomes that were collaborative samples processed alongside our pipeline (s353:76%, s354:47%, 110-FaquH:3%)
cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt


######

#!/bin/bash -l
#SBATCH -J filt_AN_97inds
#SBATCH --account=project_2001443
#SBATCH -o filt_AN_97inds.out
#SBATCH -e filt_AN_97inds.err
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --ntasks 4
#SBATCH --mem=8G

# FILTERING BASED ON ALLELIC NUMBER (AN)
# 93 diploid samples = 186 alleles expected
# 10% missing data = 167,4 = 168 alleles min.

VCF1=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.vcf.gz
VCF2=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.93inds.vcf.gz
VCF3=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.93inds.AN10percMiss.vcf.gz


#filter out individuals RN418, 105-FaquH, 54-Frufa, s353, s354, 110-FaquH, RN421, RN425, RN426, RN422
vcftools --gzvcf $VCF1 --remove-indv 110-FaquH  --remove-indv RN418 --remove-indv 105-FaquH --remove-indv 54-Frufa --remove-indv s353 --remove-indv s354 --remove-indv RN421 --remove-indv RN425 --remove-indv RN426 --remove-indv RN422 --recode --recode-INFO-all --stdout | bgzip > $VCF2

bcftools index -t $VCF2
echo "number of inds when should be 93..."
bcftools query -l $VCF2 | wc -l
echo "number of variants when should be 3.662.669..."
bcftools index -n $VCF2

#Filter on allelic number (AN)
bcftools view --threads 4 -Oz -i 'AN >= 168' $VCF2 > $VCF3

bcftools index -t $VCF3
echo "number of variants after AN filtering..."
bcftools index -n $VCF3

###END.

######################################################## Beginning script 5


#Exclude Scaffold03 (and indicate in the filename it has no Scaffold00 anyway) - for all analyses at this stage as we look for genome-wide correlations between variables

sinteractive...
cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt
module load biokit

VCFIN=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.93inds.AN10percMiss.vcf.gz

vcftools --gzvcf $VCFIN --not-chr Scaffold03 --recode --recode-INFO-all --stdout | bgzip > DP8.93inds.AN10.noScaff0003.vcf.gz       # Not directly for any analyses; 2.556.096 variants
bcftools index -t DP8.93inds.AN10.noScaff0003.vcf.gz 

### Run the rest as a batch job filter_mac_thin.sh ###

#!/bin/bash -l
#SBATCH -J filter_mac_thin
#SBATCH --account=project_2001443
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/filt/logs/filter_mac_thin.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/vcf/filt/logs/filter_mac_thin.err
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --ntasks 4
#SBATCH --mem=8G
#SBATCH --mail-type=END

cd /scratch/project_2001443/barriers_introgr_formica/vcf/filt
module load biokit

vcftools --gzvcf DP8.93inds.AN10.noScaff0003.vcf.gz --mac 2 --recode --recode-INFO-all --stdout | bgzip >  DP8.93inds.AN10.noScaff0003.mac2.vcf.gz       # Then thin w 1kb, add exsecta, select inds -> NJtree; thin w 20kb -> neighbournet
bcftools index -t DP8.93inds.AN10.noScaff0003.mac2.vcf.gz 
echo "number of variants in DP8.93inds.AN10.noScaff0003.mac2.vcf.gz..."
bcftools index -n DP8.93inds.AN10.noScaff0003.mac2.vcf.gz 

vcftools --gzvcf DP8.93inds.AN10.noScaff0003.mac2.vcf.gz --thin 20000 --recode --recode-INFO-all --stdout | bgzip >  DP8.93inds.AN10.noScaff0003.mac2.thin20kb.vcf.gz       # Neighbournet
bcftools index -t DP8.93inds.AN10.noScaff0003.mac2.thin20kb.vcf.gz
echo "number of variants in DP8.93inds.AN10.noScaff0003.mac2.thin20kb.vcf.gz..."
bcftools index -n DP8.93inds.AN10.noScaff0003.mac2.thin20kb.vcf.gz


vcftools --gzvcf DP8.93inds.AN10.noScaff0003.mac2.vcf.gz --thin 1000 --recode --recode-INFO-all --stdout | bgzip >  DP8.93inds.AN10.noScaff0003.mac2.thin1kb.vcf.gz       # NJtree (still to add exsecta and select inds)
bcftools index -t DP8.93inds.AN10.noScaff0003.mac2.thin1kb.vcf.gz  
echo "number of variants in DP8.93inds.AN10.noScaff0003.mac2.thin1kb.vcf.gz  ..."
bcftools index -n DP8.93inds.AN10.noScaff0003.mac2.thin1kb.vcf.gz  

###END.


