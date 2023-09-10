### Create a list of scaffolds

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit
for i in {01..27}; do echo "Scaffold$i"; done > scaffold.list

### Run a batch job to phase the data with ShapeIt4

#!/bin/bash -l
#SBATCH -J shpt4
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/logs/shpt4_Scaffold%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/logs/shpt4_Scaffold%a.err
#SBATCH --account=project_2001443
#SBATCH -t 06:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 4
#SBATCH --mem=8G

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

#Path for the SHAPEIT4 executable
SHPT4PATH=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/shapeit4_package/shapeit4-4.2.2/bin

mychr=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaffold.list)

$SHPT4PATH/shapeit4.2 --input ../whatshap/all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.vcf.gz \
  --output "all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.shapeit."$mychr".vcf.gz" \
  --region $mychr --sequencing --thread 4 \
  --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --use-PS 0.0001

### End of bash script
5_shapeit.sh (END)



------


### In case there's a need to update to ShapeIt5, here's a script to modify further ###

#!/bin/bash -l
#SBATCH -J shpt5
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/logs/shpt5_Scaffold%a.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/logs/shpt5_Scaffold%a.err
#SBATCH --account=project_2001443
#SBATCH -t 06:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 4
#SBATCH --mem=8G

cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit

mychr=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaffold.list)

export PATH="/projappl/project_2001443/shapeit5env/bin:$PATH" 

phase_common --input ../whatshap/all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.vcf.gz \
  --output "all_samples.DP8.hwe.AN10.noScaff00.mac2.whap.shapeit."$mychr".vcf.gz" \
  --region $mychr --sequencing --thread 4 \
  --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --use-PS 0.0001

### End of bash script
