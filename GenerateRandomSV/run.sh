###########################################
# Hg19 genome assembly used, replace k100.umap.bed if you use Hg38
###########################################

module load gcc/6.2.0 R/3.6.1 bedtools/2.29.0

#script
generate_random=/path/to/00-generate-random-sv/generate-random-sv.Rscript


# simple-sv.bedpe can be a file with all simple SVs from all samples, here only use SVs of PD10010a as an example
samplei=PD10010a

# generate 4 SVs for each real/observed SV, you can use any positive integer
Rscript $generate_random $samplei /path/to/simple-sv.bedpe /path/to/output /path/to/k100.umap.bed 4

# merge real and random sv together in file real_and_random.txt
awk 'FNR==1{if (NR==1) print  $0; next} {print $0}' $output/*.tsv > $output/real_and_random.txt

# make a breakpoint table
cut -f 1,2,3,4,8,9,10,12,13,14,15,16,17,18 $output/real_and_random.txt > $output/real_random_brkpt.txt
tail -n+2 $output/real_and_random.txt | cut -f 1,5,6,7,8,9,11,12,13,14,15,16,17,18 >> $output/real_random_brkpt.txt

# add unique "chr" order
awk -v FS='\t' -v OFS='\t' 'NR>1{$2="chr"$2;$NF=$NF"\t"NR-1}NR==1{$0=$0"\t""order"} 1' $output/real_random_brkpt.txt > $output/real_random_brkpt_withorder.txt

# create a bed file
awk -v FS='\t' -v OFS="\t" 'FNR!=1 {print $2,$3,$3,$7,$1,$5,$6,$8,$9,$10,$11,$12,$13,$14,$15}' $output/real_random_brkpt_withorder.txt > $output/real_random_brkpt_withorder.bed

# add gene name and ensg id
bedtools intersect -a $output/real_random_brkpt_withorder.bed -b  $input/coding_gene_full.bed -wa -wb -loj  >  $output/real_random_brkpt_withorder.gene.bed
