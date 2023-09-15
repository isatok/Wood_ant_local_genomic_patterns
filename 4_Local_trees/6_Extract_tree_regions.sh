#Speed up TWISST analysis: run it only for a specific scaffold or region. 
# a) save the header line in the window data file,
# b) grep specific scaffold to be used,
# c) join the header and the extracted region
# d) pick the same amount of lines from the trees file (has no header),
# e) run TWISST with this reduced faster-to-run dataset.

### MODIFY FROM HERE ON FOR THE CURRENT DATA AND ANALYSIS ###

WINDOWDATA=/scratch/project_2001443/analysis/twisst/trees/exs/phased_exs_out_max50percMi.phyml_bionj.w50.data.tsv


# a) Save the header ("chom", "start", "end") in the window data
head -1 $WINDOWDATA > windowdata_header.txt

# b) Extract the required region from the window data
cat $WINDOWDATA | grep "Scaffold01" > windowdata_scaff01.tsv

# c) Join the header and the extracted region
cat windowdata_header.txt windowdata_scaff01.tsv > phased_exs_out_max50percMi.phyml_bionj.w50.data.scaff01.tsv

#How many lines the windowdata from the extracted region has?
wc -l windowdata_scaff01.tsv #2249

# d) Extract the same number of lines from the tree file (has no header)
zcat phased_exs_out_max50percMi.phyml_bionj.w50.trees.gz | head -2249 | bgzip > phased_exs_out_max50percMi.phyml_bionj.w50.trees.scaff01.gz


#out files
phased_exs_out_max50percMi.phyml_bionj.w50.data.scaff01.tsv
phased_exs_out_max50percMi.phyml_bionj.w50.trees.scaff01.gz
#in
/scratch/project_2001443/analysis/twisst/trees/exs

###
### ------- OR combine chr-level tree results (as from Simon Martin's sticcs.py script output) into whole-genome tree file -------
###

WINDOWDATA=

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
  cat "phased_exs_out.sticcs.Scaffold10.trees.gz";
done > phased_exs_out.sticcs.trees.gz
#####HERE TAKE CARE OF THE COMPRESSION????######


#out files
phased_exs_out_max50percMi.phyml_bionj.w50.data.scaff01.tsv
phased_exs_out_max50percMi.phyml_bionj.w50.trees.scaff01.gz
#in
/scratch/project_2001443/analysis/twisst/trees/exs
