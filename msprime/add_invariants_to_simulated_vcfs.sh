### Add invariants to the msprime-simulated vcf files in order to achieve correct PIXY Dxy & Pi estimates for the simulated data. ###
### modified from Kieran Samuk's scripts at ###
### https://github.com/ksamuk/pixy_analysis/blob/main/data_generation/01_simulating-test-data/scripts/inject_invariant_sites.sh ###

### Choose the right working directory based on the simulation
simdir=/scratch/project_2001443/barriers_introgr_formica/msprime/aquFI_polFI_sim_1
#simdir=/scratch/project_2001443/barriers_introgr_formica/msprime/aquSWI_polWSWI_sim_1

vcf=output_corrHead.vcf.gz

vcfdir=$simdir
outdir=$simdir/invariants

cd $simdir

mkdir inject_tmp
mkdir $outdir

### list the vcf files
ls $vcfdir/output_corrHeadPos.vcf.gz > inject_tmp/vcf_files.tmp

### exemplar vcf file for building a fake invariants sites vcf; extract the first vcf from the list we created
vcfex=$(cat inject_tmp/vcf_files.tmp | head -n 1) 

### declare the total length of simulated data and block number and length: block_length * n(contig)
block_length=$(gunzip -c $vcfex | grep length | sed 's/.*length=//g' | sed 's/>//g' | head -1) #10 000
seq_number=$(gunzip -c $vcfex | grep length | sed 's/.*length=//g' | sed 's/>//g' | wc -l) #10
let seq_length=block_length*seq_number

### count the number of samples in the file
n_samples=$(gunzip -c $vcfex | grep "#CHROM" | awk '{print NF; exit}')
let n_samples=n_samples-9 # there are 9 non-genotype columns

### all the possible sites (containing both chr and position info)
rm inject_tmp/all_sites.tmp #remove if it exists

for ((scaffold=1; scaffold<=$seq_number; scaffold++)); do
  for ((number=1; number<=$block_length; number++)); do
    echo "Scaffold$scaffold $number" >> inject_tmp/all_sites.tmp
  done
done


###

#while read vcf
#do

### takes a msprime vcf as input. give the outfile a name
outfile=$(echo $vcf| sed 's/.vcf.gz$/_invar.vcf.gz/g')

echo "injecting invariant sites into $vcf..."
echo "will write to $outfile..."
echo "$n_samples samples"
echo "$seq_length sites"

### List sites with variants (containing chr info)
gunzip -c $vcf | grep -v "#" | awk '{print $1, $2}' > inject_tmp/var_sites.tmp 

### Find sites that need an invariant
diff --new-line-format="" --unchanged-line-format="" inject_tmp/all_sites.tmp inject_tmp/var_sites.tmp > inject_tmp/invar_sites.tmp
invar_sites=$(cat inject_tmp/invar_sites.tmp)

### Create a VCF with all invariant sites

### the start of a blank row
row=".\tA\t.\t.\tPASS\t.\tGT" #Use "A" for the REF base in invariant sites (must be A/T/C/G for PIXY to accept it) and "." for the ALT base.
rm inject_tmp/vcf_blank_spaces.vcf #remove if it exists

while read site
do
	gt=$(printf '0|0\t%.0s' $(eval echo "{1..$n_samples}"))

	newline="$site\t$row\t$gt"
	echo -e $newline >> inject_tmp/vcf_blank_spaces.vcf

done < inject_tmp/invar_sites.tmp

rm inject_tmp/invar_sites.tmp
rm inject_tmp/var_sites.tmp

### Substitute spaces with tabs and remove leading and trailing spaces
sed 's/\s/\t/g' inject_tmp/vcf_blank_spaces.vcf | sed 's/^[ \t]*//;s/[ \t]*$//' > inject_tmp/vcf_blank.vcf

gunzip -c $vcf | grep "#" > inject_tmp/vcf_header.vcf

gunzip -c $vcf | grep -v "#" > inject_tmp/vcf_variants.vcf

cat inject_tmp/vcf_blank.vcf inject_tmp/vcf_variants.vcf > inject_tmp/test_tmp.vcf

cat inject_tmp/test_tmp.vcf | sort -k1V,1 -k2n,2 > inject_tmp/test.vcf

cat inject_tmp/vcf_header.vcf inject_tmp/test.vcf | grep . | bgzip -c > $outfile 
tabix -p vcf $outfile 
rm inject_tmp/vcf_header.vcf inject_tmp/vcf_blank_spaces.vcf inject_tmp/vcf_variants.vcf inject_tmp/test_tmp.vcf inject_tmp/test.vcf

echo "wrote to $outfile"

done < inject_tmp/vcf_files.tmp

#clean up
cp $vcfdir/*_invar.vcf.gz $outdir

rm $vcfdir/*_invar.vcf.gz 

rm -r inject_tmp
