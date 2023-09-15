# This script is for 1) extracting regions from a local tree whole-genome file (useful for e.g. PhyML trees) or
# for 2) combining regions in cases when local trees are in separate per-chromosome files (e.g. sticcs trees).

###
### 1) Extract
###

#Speed up TWISST analysis: run it only for a specific scaffold or region. 
# a) Save the header in the window data
# b) Extract the required region from the window data
# c) Join the header and the extracted region
# d) Extract the same number of lines from the tree file (has no header)

cd /scratch/project_2001443/barriers_introgr_formica/local_trees/phyml/
WINDOWDATA=phased_exs_out_max50percMi.phyml_bionj.w50.data.tsv
REGION=Scaffold01 #CHANGE THE REGION TO DESIRED

# a) Save the header ("scaffold", "start", "end", "mid", "sites", "lnL") in the window data
head -1 $WINDOWDATA > windowdata_header.txt

# b) Extract the required region from the window data
cat $WINDOWDATA | grep "$REGION" > windowdata_$REGION.tsv

# c) Join the header and the extracted region
cat windowdata_header.txt windowdata_$REGION.tsv > phased_exs_out_max50percMi.phyml_bionj.w50.data.$REGION.tsv

#How many lines the windowdata from the extracted region has?
WCOUNT=$(wc -l < windowdata_$REGION.tsv)  #795

# d) Extract the same number of lines from the tree file (has no header)
zcat phased_exs_out_max50percMi.phyml_bionj.w50.trees.gz | head -$WCOUNT | bgzip > phased_exs_out_max50percMi.phyml_bionj.w50.trees.$REGION.gz

# Remove obsolete files
rm windowdata_header.txt windowdata_Scaffold01.tsv

### Results in out files
phased_exs_out_max50percMi.phyml_bionj.w50.data.$REGION.tsv
phased_exs_out_max50percMi.phyml_bionj.w50.trees.$REGION.gz

### Located in 
/scratch/project_2001443/barriers_introgr_formica/local_trees/phyml/


###
### ------- OR combine chr-level tree results (as from Simon Martin's sticcs.py script output) into whole-genome tree file -------
###

cd /scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs/
WINDOWDATA=phased_exs_out.sticcs.Scaffold01.window_data.tsv

# a) Save the header ("chom", "start", "end") in the window data of any windowdata file
head -1 $WINDOWDATA > windowdata_header.txt

# b) paste the separate windowdata files together in the scaffold order (without the headers)

for i in {01..27};
do
  tail -n +2 "phased_exs_out.sticcs.Scaffold${i}.window_data.tsv";
done > phased_exs_out.sticcs.window_data.tsv.tmp

# c) paste the separate header and the (all scaffolds) windowdata files together

cat windowdata_header.txt phased_exs_out.sticcs.window_data.tsv.tmp > phased_exs_out.sticcs.window_data.tsv
rm phased_exs_out.sticcs.window_data.tsv.tmp

# d) paste then also all tree files together in the scaffold order
for i in {01..27};
do
  cat "phased_exs_out.sticcs.Scaffold${i}.trees.gz";
done > phased_exs_out.sticcs.trees.gz


### Results in out files
phased_exs_out.sticcs.window_data.tsv
phased_exs_out.sticcs.trees.gz

### Located in 
/scratch/project_2001443/barriers_introgr_formica/local_trees/sticcs/

