#TREES
#containerize this in Puhti:
export PATH="/projappl/project_2001443/localgnm/bin:$PATH"
#it has, in addition to bioinfo_1222_env, phyml & tabix.

---

#SNP PIPELINE

conda create -n bioinfo_1222_env
conda activate bioinfo_1222_env

conda install -c conda-forge -c bioconda -c defaults vcflib mosdepth vt bamutil  

conda env export --from-history -n bioinfo_1222_env | grep -v prefix > env.yml

export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

---

#PIXY

conda create -n pixyenv
conda activate pixyenv

conda install -c conda-forge pixy

#usage:
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

---

#TWISST

conda create -n ete3env python=3.4.10
conda activate ete3env
conda install -c etetoolkit ete3 ete_toolchain
conda env export --from-history -n ete3env | grep -v prefix > env.yml

#usage:
export PATH="/projappl/project_2001443/ete3env/bin:$PATH"

---

#PHASING 1/2 WHATSHAP

conda create -n whatshapenv 
conda activate whatshapenv
conda install whatshap \
--channel conda-forge \
  --channel bioconda \
  --channel defaults \
  --strict-channel-priority

whatshap --version #should be v1.7 (2022-12-01) as of 28.04.2023

conda env export --from-history -n whatshapenv | grep -v prefix > env.yml

conda create -n myenv samtools bwa \
  --channel conda-forge \
  --channel bioconda \
  --channel defaults \
  --strict-channel-priority

---

#PHASING 2/2 SHAPEIT5

see when relevant: https://odelaneau.github.io/shapeit5/docs/installation/build_from_source/required_libraries
