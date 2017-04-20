#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

export PATH=/home/skgthab/Programs/bedtools2.22.1/bin:$PATH
export PATH=/home/skgthab/Programs/samtools-0.1.19:$PATH
export PATH=/home/skgthab/Programs/custom_scripts:$PATH
export PATH=/home/skgthab/Programs:$PATH
export PATH=/home/skgthab/Programs/HTSeq-0.6.1/build/scripts-2.7:$PATH
export PATH=/home/skgthab/Programs/weblogo-3.3:$PATH
export PATH=/home/skgthab/Programs/Homer-v4.5/bin:$PATH
export PATH=/home/skgthab/Programs/weblogo.2.8.2:$PATH
export PATH=/home/skgthab/Programs/BEDOPSv2.4.2:$PATH
export PATH=/home/skgthab/Programs/fastx_toolkit_0.0.13:$PATH
export PATH=/home/skgthab/Programs/cufflinks-2.2.0.Linux_x86_64:$PATH
export PATH=/home/skgthab/Programs/ribopicker-standalone-0.4.3:$PATH
export PATH=/home/skgthab/Programs/sratoolkit.2.5.0-1-ubuntu64/bin:$PATH
export PATH=/home/skgthab/Programs/iCLIPro-0.1.1/scripts:$PATH

data=$1
path=`pwd -P`
path=${path}/
genome=/media/skgthab/storage/UCL/genome-data/ucsc.hg19.fasta

gunzip ${path}${data}.gz

################
## Formating ###
################

# get cross link positions and convert it to bedGraph
python ./scripts/BEDtoXlink-CLIP.py ${path}${data} ${path}${data}-xlink.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-xlink.bed > ${path}${data}-xlink-sorted.bed
python ./scripts/BEDsum.py ${path}${data}-xlink-sorted.bed ${path}${data}-xlink-sorted-sum.bed
python ./scripts/BED2BEDgraph.py ${path}${data}-xlink-sorted-sum.bed ${path}${data}-xlink-sorted-sum.bedGraph
rm ${path}${data}-xlink.bed

# print unique number of cDNAs
echo "Unique number of cDNAs" > ${path}${data}-01.Mapping-REPORT.txt
wc -l ${path}${data} >> ${path}${data}-01.Mapping-REPORT.txt

# count cDNAs that have more then 1 crosslink (cDNA start > 1)
python ./scripts/count_cDNAs.py ${path}${data}-xlink-sorted-sum.bed >> ${path}${data}-01.Mapping-REPORT.txt

# compress the original data
gzip ${path}${data}
gzip ${path}${data}-xlink-sorted-sum.bed
gzip ${path}${data}-xlink-sorted-sum.bedGraph

#################################
### kmers at crosslinks ###
#################################

# flank end positions by 2nt upostream and downstream
python ./scripts/flankBEDpositionsCustom.py ${path}${data}-xlink-sorted.bed ${path}${data}-xlink-sorted-flanked2.bed 2 2

# get fasta around cross links
bedtools getfasta -s -fi ${genome} -bed ${path}${data}-xlink-sorted-flanked2.bed -fo ${path}${data}-xlink-sorted-flanked2.fasta

# find kmers around cross links
python ./scripts/kmer_finder.py ${path}${data}-xlink-sorted-flanked2.fasta ${path}${data}-xlinks-flanked2-kmers.txt

gzip ${path}${data}-xlink-sorted.bed
rm ${path}${data}-xlink-sorted-flanked2.bed
rm ${path}${data}-xlink-sorted-flanked2.fasta

#################################
### CLIPo length distribution ###
#################################

# how narrow is the length distribution
python ./scripts/narrow_cDNA_length_distribution.py ${path}${data} >> ${path}${data}-01.Mapping-REPORT.txt
gzip ${path}${data}.bed
