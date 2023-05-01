#parceVCF.py script:
/scratch/project_2001443/analysis/genomics_simon/genomics_general/VCF_processing/parseVCF.py

#Phased vcf, with F. exsecta added in:
FULLVCF=/scratch/project_2001443/vcf/phasing/shapeit/phased_with_outgroup_gtfix.vcf.gz

#Sample table (still without the outgroup, F. exsecta)
/scratch/project_2001443/vcf/ind_lists/sample_table.tab

#Create sample table with F. exsecta
nano fexs.tab #fexs.tab contains three columns: "Fexs  exsecta exsecta"
cat /scratch/project_2001443/vcf/ind_lists/sample_table.tab fexs.tab > sample_table_fexs.tab

#create a ploidy table
cut -f1 sample_table_fexs.tab | grep -v 'vcfID' | sed 's/$/ 2/' > ploidy.tab
nano ploidy.tab #Change manually F. exsecta to "1"
PLOIDYTAB=/scratch/project_2001443/vcf/ind_lists/ploidy.tab

#make a geno file
cd /scratch/project_2001443/vcf/geno
PARSEPY=/scratch/project_2001443/analysis/genomics_simon/genomics_general/VCF_processing/parseVCF.py
python $PARSEPY -i $FULLVCF --skipIndels --excludeDuplicates --ploidyFile $PLOIDYTAB | bgzip > phased.exsecta.geno.gz
