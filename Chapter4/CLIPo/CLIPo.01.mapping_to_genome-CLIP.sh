#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

# set correct paths for the following tools
export PATH=~/Programs/bedtools2.22.1/bin:$PATH
export PATH=~/Programs/samtools-0.1.19:$PATH
export PATH=~/Programs/scripts:$PATH
export PATH=~/Programs/fastx_toolkit_0.0.13:$PATH

data=$1
path=`pwd -P`
path=${path}/
genome=ucsc.hg19.fasta

gunzip ${path}${data}.gz

############################
## Mapping and formating ###
############################

# clip the adapter and discard clipped sequences and discard the sequences that are shorter then 17 nt + 5 random barcode + 4 experimental barcode; complete cDNAs contains the adapter sequence and incomplete does not
fastx_clipper -Q 33 -a AGATCGGAAG -C -n -l 26 -i  ${path}${data} -o  ${path}${data}-incomplete.fq
fastx_clipper -Q 33 -a AGATCGGAAG -c -n -l 26 -i  ${path}${data} -o  ${path}${data}-complete.fq

# fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-incomplete.fq -o ${path}${data}-incomplete.fa
fastq_to_fasta -Q 33 -n -i ${path}${data}-complete.fq -o ${path}${data}-complete.fa
rm ${path}${data}-incomplete.fq
rm ${path}${data}-complete.fq

# swap random barcodes to headers of fasta file
python ./scripts/${path}swap_barcode_to_header.py ${path}${data}-incomplete.fa ${path}${data}-incomplete-barcodes.fa
python ./scripts/${path}swap_barcode_to_header.py ${path}${data}-complete.fa ${path}${data}-complete-barcodes.fa
rm ${path}${data}-incomplete.fa
rm ${path}${data}-complete.fa

# map to hg19
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-incomplete-barcodes.fa -S ${path}${data}-incomplete.sam
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-complete-barcodes.fa -S ${path}${data}-complete.sam
rm ${path}${data}-incomplete-barcodes.fa
rm ${path}${data}-complete-barcodes.fa

# filter reads with more then 2 mismatches
samtools view -Sh ${path}${data}-incomplete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-incomplete-2mis.sam
samtools view -Sh ${path}${data}-complete.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-complete-2mis.sam

# SAM to BAM
samtools view -hSb ${path}${data}-incomplete-2mis.sam > ${path}${data}-incomplete-2mis.bam
samtools view -hSb ${path}${data}-complete-2mis.sam > ${path}${data}-complete-2mis.bam
rm ${path}${data}-incomplete-2mis.sam
rm ${path}${data}-complete-2mis.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-incomplete-2mis.bam > ${path}${data}-incomplete-2mis.bed
bedtools bamtobed -i ${path}${data}-complete-2mis.bam > ${path}${data}-complete-2mis.bed

# remove duplicates
cat ${path}${data}-incomplete-2mis.bed | sort -k1,1 -k2,2n -k5,5 | > ${path}${data}-incomplete-2mis.bed
cat ${path}${data}-complete-2mis.bed | sort -k1,1 -k2,2n -k5,5 | > ${path}${data}-complete-2mis.bed
rm ${path}${data}-incomplete-2mis.bed
rm ${path}${data}-complete-2mis.bed

# merge complete and incompelte cDNAs
cat ${path}${data}-incomplete-2mis.bed ${path}${data}-complete-2mis.bed | sort -k1,1 -k2,2n -k6,6 > ${path}${data}-all-2mis.bed
gzip ${path}${data}-incomplete-2mis.bed
gzip ${path}${data}-complete-2mis.bed

# get cross link positions and convert it to bedGraph
python ./scripts/BEDtoXlink-CLIP.py ${path}${data}-all-2mis.bed ${path}${data}-all-2mis-xlink.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-all-2mis-xlink.bed > ${path}${data}-all-2mis-xlink-sorted.bed
python ./scripts/BEDsum.py ${path}${data}-all-2mis-xlink-sorted.bed ${path}${data}-all-2mis-xlink-sorted-sum.bed
python ./scripts/BED2BEDgraph.py ${path}${data}-all-2mis-xlink-sorted-sum.bed ${path}${data}-all-2mis-xlink-sorted-sum.bedGraph
rm ${path}${data}-all-2mis-xlink.bed

# print unique number of cDNAs
echo "Unique number of cDNAs" > ${path}${data}-01.Mapping-REPORT.txt
wc -l ${path}${data}-all-2mis.bed >> ${path}${data}-01.Mapping-REPORT.txt

# count cDNAs that have more then 1 crosslink (cDNA start > 1)
python ./scripts/count_cDNAs.py ${path}${data}-all-2mis-xlink-sorted-sum.bed >> ${path}${data}-01.Mapping-REPORT.txt

# compress the original data
gzip ${path}${data}
gzip ${path}${data}-all-2mis-xlink-sorted-sum.bed
gzip ${path}${data}-all-2mis-xlink-sorted-sum.bedGraph

#################################
### kmers at crosslinks ###
#################################

# flank end positions by 2nt upostream and downstream
python ./scripts/flankBEDpositionsCustom.py ${path}${data}-all-2mis-xlink-sorted.bed ${path}${data}-all-2mis-xlink-sorted-flanked2.bed 2 2

# get fasta around cross links
bedtools getfasta -s -fi ${genome} -bed ${path}${data}-all-2mis-xlink-sorted-flanked2.bed -fo ${path}${data}-all-2mis-xlink-sorted-flanked2.fasta

# find kmers around cross links
python ./scripts/kmer_finder.py ${path}${data}-all-2mis-xlink-sorted-flanked2.fasta ${path}${data}-xlinks-flanked2-kmers.txt

gzip ${path}${data}-all-2mis-xlink-sorted.bed
rm ${path}${data}-all-2mis-xlink-sorted-flanked2.bed
rm ${path}${data}-all-2mis-xlink-sorted-flanked2.fasta

#################################
### CLIPo length distribution ###
#################################

# how narrow is the length distribution
echo "B:  Length constraint (% of cDNAs in 10nt length window)" > ${path}${data}-01.Mapping-REPORT.txt
python ./scripts/narrow_cDNA_length_distribution.py ${path}${data}-all-2mis.bed >> ${path}${data}-01.Mapping-REPORT.txt
gzip ${path}${data}-all-2mis.bed
