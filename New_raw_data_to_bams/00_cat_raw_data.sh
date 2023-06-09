#Here the purpose is to merge raw data so that each sample has only one raw data fastq file. For some samples
#there are originally two, because they needed to be sequenced more.

B1PATH=/scratch/project_2001443/barriers_introgr_formica/fastq/raw/batch1
B2PATH=/scratch/project_2001443/barriers_introgr_formica/fastq/raw/batch2
COMBPATH=/scratch/project_2001443/barriers_introgr_formica/fastq/raw/combined

#batch2 contains two files per sample (L2 & L3), and the L3 is identical to the file delivered already in batch1
#checked this with md5sums, e.g., 
md5sum $B1PATH/RN416/*_1.fq.gz $B2PATH/RN416/*_L3_1.fq.gz
#e9156ccbc1121a35f3333d508a6606bc  /scratch/project_2001443/barriers_introgr_formica/fastq/raw/batch1/RN416/RN416_1.fq.gz
#e9156ccbc1121a35f3333d508a6606bc  /scratch/project_2001443/barriers_introgr_formica/fastq/raw/batch2/RN416/RN416_EKDL230000676-1A_HNHM2DSX5_L3_1.fq.gz

#combine for samples 415,416,417 two different files provided in batch2, else use the one in batch1
cat $B2PATH/RN415/*_L3_1.fq.gz $B2PATH/RN415/*_L2_1.fq.gz > $COMBPATH/RN415_comb_1.fq.gz
cat $B2PATH/RN416/*_L3_1.fq.gz $B2PATH/RN416/*_L2_1.fq.gz > $COMBPATH/RN416_comb_1.fq.gz
cat $B2PATH/RN417/*_L3_1.fq.gz $B2PATH/RN417/*_L2_1.fq.gz > $COMBPATH/RN417_comb_1.fq.gz

cat $B2PATH/RN415/*_L3_2.fq.gz $B2PATH/RN415/*_L2_2.fq.gz > $COMBPATH/RN415_comb_2.fq.gz
cat $B2PATH/RN416/*_L3_2.fq.gz $B2PATH/RN416/*_L2_2.fq.gz > $COMBPATH/RN416_comb_2.fq.gz
cat $B2PATH/RN417/*_L3_2.fq.gz $B2PATH/RN417/*_L2_2.fq.gz > $COMBPATH/RN417_comb_2.fq.gz

cp $B1PATH/RN418/*fq.gz $COMBPATH/
cp $B1PATH/RN419/*fq.gz $COMBPATH/
cp $B1PATH/RN420/*fq.gz $COMBPATH/
cp $B1PATH/RN421/*fq.gz $COMBPATH/
cp $B1PATH/RN422/*fq.gz $COMBPATH/
cp $B1PATH/RN423/*fq.gz $COMBPATH/
cp $B1PATH/RN424/*fq.gz $COMBPATH/
cp $B1PATH/RN425/*fq.gz $COMBPATH/
cp $B1PATH/RN426/*fq.gz $COMBPATH/
