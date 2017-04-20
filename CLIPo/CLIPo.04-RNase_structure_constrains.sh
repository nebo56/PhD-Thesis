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
export PATH=/home/skgthab/programs/RNAfold/usr/bin:$PATH

cDNAs=$1	#complete cDNAs only
path=/SAN/neuroscience/TE/Nejc/CLIPo-all/
genome=/home/skgthab/bowtie-indexes/hg19/hg19.fa

gunzip ${path}${cDNAs}.gz

# get cDNA end
python ${path}scripts/getEnd-BED.py ${path}${cDNAs} tmp.cDNA.end.bed

# flank 30 nt in both directions
python ${path}scripts/flankBEDpositionsCustom.py tmp.cDNA.end.bed tmp.cDNA.end.flanked30.bed 30 30 
sort -k1,1 -k2,2n -k6,6 tmp.cDNA.end.flanked30.bed > tmp.cDNA.end.flanked30-sorted.bed

# get fasta
bedtools getfasta -s -fi ${genome} -bed tmp.cDNA.end.flanked30-sorted.bed -fo tmp.cDNA.end.flanked30.fasta
rm tmp.cDNA.end.flanked30-sorted.bed

# get structure from RNAfold
RNAfold --noPS -i tmp.cDNA.end.flanked30.fasta > tmp.cDNA.end.flanked30.RNAfold.fasta

# sum of pairing probability 
echo "Structure constraints at cDNA-ends" > ${path}${cDNAs}-04.CLIPo-cDNAend-structure-REPORT.txt
python ${path}scripts/sumRNAfold.py tmp.cDNA.end.flanked30.RNAfold.fasta 27 32 >> ${path}${cDNAs}-04.CLIPo-cDNAend-structure-REPORT.txt

#clean
rm tmp*

########### control ########### 

# get cDNA end
python ${path}scripts/getEnd-BED.py ${path}${cDNAs} tmp.cDNA.end.bed

# get -10 nt region from cDNA end position
python ${path}scripts/flankBEDpositionsCustom.py tmp.cDNA.end.bed tmp.cDNA.end-10.bed 10 -10

# flank 50 nt in both directions
python ${path}scripts/flankBEDpositionsCustom.py tmp.cDNA.end-10.bed tmp.cDNA.end.flanked30.bed 30 30
sort -k1,1 -k2,2n -k6,6 tmp.cDNA.end.flanked30.bed > tmp.cDNA.end.flanked30-sorted.bed

# get fasta
bedtools getfasta -s -fi ${genome} -bed tmp.cDNA.end.flanked30-sorted.bed -fo tmp.cDNA.end.flanked30.fasta
rm tmp.cDNA.end.flanked30-sorted.bed

# get structure from RNAfold
RNAfold --noPS -i tmp.cDNA.end.flanked30.fasta > tmp.cDNA.end.flanked30.RNAfold.fasta

# sum of pairing probability 
echo "Structure constraints at cDNA-ends (control -10nt)" >> ${path}${cDNAs}-04.CLIPo-cDNAend-structure-REPORT.txt
python ${path}scripts/sumRNAfold.py tmp.cDNA.end.flanked30.RNAfold.fasta 27 32 >> ${path}${cDNAs}-04.CLIPo-cDNAend-structure-REPORT.txt

#clean
rm tmp*
