#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

export PATH=~/Programs/bedtools2.22.1/bin:$PATH
export PATH=~/Programs/scripts:$PATH
export PATH=~/Programs/RNAfold/usr/bin:$PATH

cDNAs=$1	
path=`pwd -P`
path=${path}/
genome=hg19.fasta

# separate cDNAs that are shorter than 40 nt
python ./scripts/separate_cDNAs_40nt_length.py ${path}${data} ${path}${data}-40less.bed

# get cDNA end
python ${path}scripts/getEnd-BED.py ${path}${cDNAs}-40less.bed tmp.cDNA.end.bed

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

# clean up
rm tmp*
