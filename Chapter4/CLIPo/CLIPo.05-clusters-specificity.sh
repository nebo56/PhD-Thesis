#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

export PATH=~/Programs/bedtools2.22.1/bin:$PATH
export PATH=~/Programs/scripts:$PATH
export PATH=~/Programs/RNAfold/usr/bin:$PATH

data=$1
path=`pwd -P`
path=${path}/
genome=hg19.fasta
clusters=$2


### find kmers ###

# get start and end positions
python ${path}scripts/getStartAndEnd-BED.py ${path}${data} ${path}${data}

# flank end position by 5 nt upostream and downstream
python ${path}scripts/flankBEDpositionsCustom.py ${path}${data}-Start.bed ${path}${data}-Start-flanked5.bed 5 5

# get fasta around cDNA start positions
bedtools getfasta -s -fi ${genome} -bed ${path}${data}-Start-flanked5.bed -fo ${path}${data}-Start-flanked5.fasta

# find kmers around cDNA starts
python ${path}scripts/kmer_finder.py ${path}${data}-Start-flanked5.fasta ${path}${data}-Start-flanked5-kmers.txt



### kmers cluster enrichment ###

# report
echo "Number of binding sites (clusters)" > ${path}${data}-03.REPORT-kmer_enrichment.txt
wc -l ${path}${clusters} >> ${path}${data}-03.REPORT-kmer_enrichment.txt

#ignore clusters that are smaller then 5nt
python ${path}scripts/remove_clusters_5less.py ${path}${clusters} tmp.clusters.5more.bed

# get fasta
bedtools getfasta -s -fi ${genome} -bed tmp.clusters.5more.bed -fo tmp.clusters.5more.fasta

echo "Top 10 kmers from (-5..X..+5) coverage inside of clusters" >> ${path}${data}-03.REPORT-kmer_enrichment.txt
python ${path}scripts/kmer_coverage-top10.py tmp.clusters.5more.fasta ${path}${data}-Start-flanked5-kmers.txt >> ${path}${data}-03.REPORT-kmer_enrichment.txt

# get control outside from cluster region
python ${path}scripts/getStartAndEnd-BED.py ${path}${clusters} tmp.clusters
python ${path}scripts/flankBEDpositionsCustom.py tmp.clusters-End.bed tmp.clusters-End-flanked-500_600.bed -500 600

# remove controls that are overlapping clusters
bedtools intersect -s -v -a tmp.clusters-End-flanked-500_600.bed -b tmp.clusters.5more.bed > tmp.clusters-End-flanked-500_600-no_overlaps.bed

# get fasta
bedtools getfasta -s -fi ${genome} -bed tmp.clusters-End-flanked-500_600-no_overlaps.bed -fo tmp.clusters-End-flanked-500_600-no_overlaps.fasta

# report
echo "- Kmer enrichment in clusters (control -500..cDNAend..-600)" >> ${path}${data}-03.REPORT-kmer_enrichment.txt
python ${path}scripts/kmer_coverage-top10.py tmp.clusters-End-flanked-500_600-no_overlaps.fasta ${path}${data}-Start-flanked5-kmers.txt >> ${path}${data}-03.REPORT-kmer_enrichment.txt

#clean up
rm tmp
gzip ${path}${data}
gzip ${path}${clusters}

