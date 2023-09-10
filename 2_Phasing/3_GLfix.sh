#!/bin/bash -l
#SBATCH -J GLfix
#SBATCH -o /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap/logs/GLfix.out
#SBATCH -e /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap/logs/GLfix.err
#SBATCH --account=project_2001443
#SBATCH -t 08:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=6G
#SBATCH --mail-type=END


cd /scratch/project_2001443/barriers_introgr_formica/vcf/phasing/whatshap
module load biokit

for i in *.phased.vcf.gz ; do
  echo $i
  mv $i ${i}.tmp
  bcftools view ${i}.tmp -Ov | perl -pi -e 's/(\#\#FORMAT=<ID=GL,Number=)G/$1\./' | bgzip > $i
  bcftools index -t $i
done

###end
3_GLfix.sh (END)
