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
#path=`pwd -P`
path=/SAN/neuroscience/TE/Nejc/CLIPo-all/
genome=/home/skgthab/bowtie-indexes/hg19/hg19.fa

gunzip ${path}${data}.gz	#you should add complete cDNAs 

################
### CLIPo 40 ###
################

# get cDNAs shorter then 40 from complete ones
python ${path}/scripts/separate_cDNAs_40nt_length.py ${path}${data} ${path}${data}-40less.bed
echo "A: number of uniquely mapped cDNAs " > ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt
wc -l ${path}${data} >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt
echo "A: number of uniquely mapped cDNAs (< 40nt length)" >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt
wc -l ${path}${data}-40less.bed >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt

gzip ${path}${data} 
data=${data}-40less.bed

# get cDNAs that have dominant start, dominant end or equal
#bash separate-cDNAs.sh ${path}${data}-40less.bed
#echo "Dominant start cDNAs (<40) " >> ${path}${data}-02.CLIPo-REPORT.txt
#wc -l ${path}${data}-40less.bed-cDNAstartScore.bed >> ${path}${data}-02.CLIPo-REPORT.txt
#echo "Dominant end cDNAs (<40) " >> ${path}${data}-02.CLIPo-REPORT.txt
#wc -l ${path}${data}-40less.bed-cDNAendScore.bed >> ${path}${data}-02.CLIPo-REPORT.txt
#echo "Equal cDNAs (<40) " >> ${path}${data}-02.CLIPo-REPORT.txt
#wc -l ${path}${data}-40less.bed-cDNAequalScore.bed >> ${path}${data}-02.CLIPo-REPORT.txt

#gzip ${path}${data}-40less.bed-cDNAscore.bed
#gzip ${path}${data}-40less.bed-cDNAstartScore.bed
#gzip ${path}${data}-40less.bed-cDNAendScore.bed
#gzip ${path}${data}-40less.bed-cDNAequalScore.bed
#gzip ${path}${data}-40less.bed

###################################
### CLIPo kmers end constraints ###
###################################

# get start and end positions
python ${path}scripts/getStartAndEnd-BED.py ${path}${data} ${path}${data}

# flank end positions by 2nt upostream and downstream
python ${path}scripts/flankBEDpositionsCustom.py ${path}${data}-End.bed ${path}${data}-End-flanked2.bed 2 2

# get fasta around cDNA ends
bedtools getfasta -s -fi ${genome} -bed ${path}${data}-End-flanked2.bed -fo ${path}${data}-End-flanked2.fasta

# find kmers around cDNA ends
python ${path}scripts/kmer_finder.py ${path}${data}-End-flanked2.fasta ${path}${data}-End-kmers.txt

# count how many cDNAs contains top 10 kmers at the end
echo "C: Top 10 kmers at cDNAend (cDNAs <40 nt)" >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt
python ${path}scripts/kmer-separation-counter.py ${path}${data}-End-flanked2.fasta ${path}${data}-End-kmers.txt >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt

# get a control fasta from region -15 to -5 upstream from cDNAend and count how mand cDNAs conatins top 10 kmers in the control region
python ${path}scripts/flankBEDpositionsCustom.py ${path}${data}-End.bed ${path}${data}-End-flanked15-5.bed 15 -5
bedtools getfasta -s -fi ${genome} -bed ${path}${data}-End-flanked15-5.bed -fo ${path}${data}-End-flanked15-5.fasta
echo "C: Contains top 10 kmers at control region -15..-5 upstream from cDNAend (cDNAs <40 nt)" >> ${path}${data}-02.CLIPokmer-constrains--REPORT.txt
python ${path}scripts/kmer-separation-counter.py ${path}${data}-End-flanked15-5.fasta ${path}${data}-End-kmers.txt >> ${path}${data}-02.CLIPokmer-constrains--REPORT.txt

# cDNA length constraints; what is the maximum distribution od cDNA lenghts in 10 nt window
echo "D: cDNA end Contains; distribution od cDNA lenghts in 10 nt window" >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt
python ${path}scripts/CLIPo_cDNA_length_distribution.py ${path}${data} >> ${path}${data}-02.CLIPo-kmer-constrains-REPORT.txt

# 

rm ${path}${data}-End-flanked2.fasta
rm ${path}${data}-End-flanked2.bed
rm ${path}${data}-End-flanked15-5.fasta
rm ${path}${data}-End-flanked15-5.bed

gzip ${path}${data}
gzip ${path}${data}-End.bed
gzip ${path}${data}-Start.bed


