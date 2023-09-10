#!/bin/bash -l
#SBATCH -J shpt4
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/logs/shpt4_Scaffold%a.out
#SBATCH -e /scratch/project_2001443/vcf/phasing/shapeit/logs/shpt4_Scaffold%a.err
#SBATCH --account=project_2001443
#SBATCH -t 06:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 4
#SBATCH --mem=8G

cd $SCRATCH/vcf/phasing/shapeit

mychr=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaffold.list)

module load bioconda/3
conda activate my_seqdata

shapeit4 --input ../whatshap/all_samples.minDP8.AN10percMiss.mac2.whap.vcf.gz \
  --output "all_samples.minDP8.AN10percMiss.mac2.whap.shapeit."$mychr".vcf.gz" \
  --region $mychr --sequencing --thread 4 \
  --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --use-PS 0.0001

### End of bash script
5_shapeit.sh (END)
