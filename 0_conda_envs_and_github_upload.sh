
### LOCAL TREES ###

#usage:
PATH="/projappl/project_2001443/phymlenv/bin:$PATH" 

## env.yml
# name: phymlenv
# channels:
#  - conda-forge
#   - bioconda
#   - defaults
# dependencies:
#   - phyml
#   - python
#   - numpy


###SEE IF THIS LOCALGNM IS NEEEDED. AT LEAST PHYML CANNOT BE RUN FROM HERE IF PYTHON/NUMPY NOT INSTALLED###
#usage:
#export PATH="/projappl/project_2001443/localgnm/bin:$PATH" 

## env.yml
# name: localgnm
# channels:
#  - conda-forge
#   - bioconda
#   - defaults
# dependencies:
#   - vt
#   - tabix
#   - bamutil
#   - mosdepth
#   - vcflib
#   - phyml

# ---------------------

### SNP PIPELINE ###

conda create -n bioinfo_1222_env
conda activate bioinfo_1222_env

conda install -c conda-forge -c bioconda -c defaults vcflib mosdepth vt bamutil  

conda env export --from-history -n bioinfo_1222_env | grep -v prefix > env.yml

#usage:
export PATH="/projappl/project_2001443/bioinfo_1222_env/bin:$PATH"

# ---------------------

### PIXY ###

conda create -n pixyenv
conda activate pixyenv

conda install -c conda-forge pixy

#usage:
export PATH="/projappl/project_2001443/pixyenv/bin:$PATH"

# ---------------------

### TWISST ###

conda create -n ete3env python=3.4.10
conda activate ete3env
conda install -c etetoolkit ete3 ete_toolchain
conda env export --from-history -n ete3env | grep -v prefix > env.yml

#usage:
export PATH="/projappl/project_2001443/ete3env/bin:$PATH"

# ---------------------

### PHASING 1/2 WHATSHAP ###

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

## env.yml
# name: whatshapenv
# channels:
#   - conda-forge
#   - bioconda
#   - defaults
# dependencies:
#   - whatshap

# ---------------------

### PHASING 2/2 SHAPEIT ###

## Using ShapeIt4 for the current analyses ##

## Install ShapeIt4 dependencies htslib and boost (suitable older versions)

#htslib version 1.9

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make

#boost version 1.67

wget https://boostorg.jfrog.io/artifactory/main/release/1.67.0/source/boost_1_67_0.tar.bz2
tar -vxjf boost_1_67_0.tar.bz2
cd boost_1_67_0
./bootstrap.sh --help
./bootstrap.sh --prefix=/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/shapeit4_package/boost_1_67_0 \
  --with-libraries=program_options,iostreams
./b2 install

## Install ShapeIt4

wget https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz
tar -vxjf v4.2.2.tar.gz
#modify requested paths in makefile to refer to the right files in packages above, and type "make"
#installed in
/scratch/project_2001443/barriers_introgr_formica/vcf/phasing/shapeit/shapeit4_package/shapeit4-4.2.2
#install shapeit5

conda create -n shapeit5env
conda activate shapeit5env
conda install -c bioconda shapeit5


## Install ShapeIt5 (version 5.1.1) - NOT USED IN THE CURRENT ANALYSES

#conda env export --from-history -n shapeit5env | grep -v prefix > env.yml

# name: shapeit5env
# channels:
#   - conda-forge
#   - bioconda
#   - defaults


#usage:
export PATH="/projappl/project_2001443/shapeit5env/bin:$PATH" 

# ---------------------

### GENOMICS_GENERAL GitHub update (last uploaded 29.06.2023) ###

cd /scratch/project_2001443/analysis/genomics_simon
git clone https://github.com/simonhmartin/genomics_general.git
