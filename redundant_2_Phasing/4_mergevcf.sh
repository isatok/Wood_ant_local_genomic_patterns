### Create a list of the phased vcf files
cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap
ls *.phased.vcf.gz > phased_vcf.list

### Run a batch job to merge the individual whatshap-phased vcf files into one vcf file

#!/bin/bash -l
#SBATCH -J mergevcf
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap/logs/mergevcf_%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap/logs/mergevcf_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 02:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=5G

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap
module load biokit

bcftools merge -l /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap/phased_vcf.list -Oz > all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.vcf.gz &&\
bcftools index -t all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.vcf.gz

###END
4_mergevcf.sh (END)
