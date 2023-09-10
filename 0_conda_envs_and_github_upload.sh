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


#8.9.2023 update bioconda set-up as instructed on the whatshap website
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

#install whatshap
conda create -n whatshapenv whatshap
conda activate whatshapenv

whatshap --version #2.0, on 08.09.2023

conda env export --from-history -n whatshapenv | grep -v prefix > env.yml

#usage:
export PATH="/projappl/project_2001443/whatshapenv/bin:$PATH" 

# name: whatshapenv
# channels:
#   - conda-forge
#   - bioconda
#   - defaults
# dependencies:
#   - whatshap


#WHY IS THIS HERE?
# conda create -n myenv samtools bwa \
#   --channel conda-forge \
#   --channel bioconda \
#   --channel defaults \
#   --strict-channel-priority

---

#PHASING 2/2 SHAPEIT


#install shapeit4

#conda create -n shapeit4env
#conda activate shapeit4env
#conda install -c bioconda shapeit4

#install shapeit5

conda create -n shapeit5env
conda activate shapeit5env
conda install -c bioconda shapeit5


#shapeit5 version 5.1.1

conda env export --from-history -n shapeit5env | grep -v prefix > env.yml

# name: shapeit5env
# channels:
#   - conda-forge
#   - bioconda
#   - defaults


#usage:
export PATH="/projappl/project_2001443/shapeit5env/bin:$PATH" 

---

#GENOMICS_GENERAL GitHub update (last uploaded 29.06.2023)

cd /scratch/project_2001443/analysis/genomics_simon
git clone https://github.com/simonhmartin/genomics_general.git
