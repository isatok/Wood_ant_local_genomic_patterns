#!/bin/bash -l
#SBATCH -J mergevcf
#SBATCH -o /scratch/project_2001443/vcf/phasing/whatshap/logs/mergevcf_%a.out
#SBATCH -e /scratch/project_2001443/vcf/phasing/whatshap/logs/mergevcf_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 02:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=5G

cd /scratch/project_2001443/vcf/phasing/whatshap
module load biokit

bcftools merge -l /scratch/project_2001443/vcf/phasing/whatshap/phased_vcf.list -Oz > all_samples.minDP8.AN10percMiss.mac2.whap.vcf.gz &&\
bcftools index -t all_samples.minDP8.AN10percMiss.mac2.whap.vcf.gz

###END
4_mergevcf.sh (END)
