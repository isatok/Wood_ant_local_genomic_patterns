#Speed up TWISST analysis: run it only for a specific scaffold or region. 
# a) save the header line in the tree windows data file,
# b) grep specific scaffold to be used,
# c) pick the same amount of lines from the trees file (has no header),
# d) run TWISST with this reduced faster-to-run dataset.

### MODIFY FROM HERE ON FOR THE CURRENT DATA AND ANALYSIS ###

TREEFILE=/scratch/project_2001443/analysis/twisst/trees/exs/phased_exs_out_max50percMi.phyml_bionj.w50.data.tsv

head -1 $TREEFILE > trees_header.txt
cat $TREEFILE | grep "Scaffold01" > trees_scaff01.tsv
cat trees_header.txt trees_scaff01.tsv > phased_exs_out_max50percMi.phyml_bionj.w50.data.scaff01.tsv

wc -l trees_scaff01.tsv #2249

zcat phased_exs_out_max50percMi.phyml_bionj.w50.trees.gz | head -2249 | bgzip > phased_exs_out_max50percMi.phyml_bionj.w50.trees.scaff01.gz


#out files
phased_exs_out_max50percMi.phyml_bionj.w50.data.scaff01.tsv
phased_exs_out_max50percMi.phyml_bionj.w50.trees.scaff01.gz
#in
/scratch/project_2001443/analysis/twisst/trees/exs
