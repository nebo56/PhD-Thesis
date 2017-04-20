#!/bin/bash -l
#
# bedtools needs to be installed to run this pipeline http://bedtools.readthedocs.io/en/latest/#
#
### input ###
clusters=$1	#crosslink clusters
exons=$2	#regulated exons
controls=$3	#control exons
kmers=$4	#motifs

### GENOME FASTA ###
fasta=ucsc.hg19.fasta	#set path to genome fasta

# get splice sites for control and regulated exons
python get3SS.py ${exons} ${exons}
python get5SS.py ${exons} ${exons}
python get3SS.py ${controls} ${controls}
python get5SS.py ${controls} ${controls}

# flank splice site region
python flankBEDpositionsCustom.py ${exons}-3SS.bed ${exons}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${exons}-5SS.bed ${exons}-5SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${controls}-3SS.bed ${controls}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${controls}-5SS.bed ${controls}-5SS-flanked300.bed 300 300

# intersect flanked exon region with clusters (includes all flanked regions)
bedtools intersect -s -b ${clusters} -a ${exons}-3SS-flanked300.bed -wao > ${exons}-3SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${exons}-5SS-flanked300.bed -wao > ${exons}-5SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${controls}-3SS-flanked300.bed -wao > ${controls}-3SS-flanked300-clusters.bed
bedtools intersect -s -b ${clusters} -a ${controls}-5SS-flanked300.bed -wao > ${controls}-5SS-flanked300-clusters.bed

# get positions of clusters relative to exons
python get_flanked_cluster-positions.py ${exons}-3SS-flanked300-clusters.bed ${exons}-3SS-flanked300-clusters-map_positions.bed
python get_flanked_cluster-positions.py ${exons}-5SS-flanked300-clusters.bed ${exons}-5SS-flanked300-clusters-map_positions.bed
python get_flanked_cluster-positions.py ${controls}-3SS-flanked300-clusters.bed ${controls}-3SS-flanked300-clusters-map_positions.bed
python get_flanked_cluster-positions.py ${controls}-5SS-flanked300-clusters.bed ${controls}-5SS-flanked300-clusters-map_positions.bed

# get fasta files of flanking exon region with cluster positions
bedtools getfasta -s -name -fi ${fasta} -bed ${exons}-3SS-flanked300-clusters-map_positions.bed -fo ${exons}-3SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${exons}-5SS-flanked300-clusters-map_positions.bed -fo ${exons}-5SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${controls}-3SS-flanked300-clusters-map_positions.bed -fo ${controls}-3SS-flanked300-clusters-map_positions.fasta
bedtools getfasta -s -name -fi ${fasta} -bed ${controls}-5SS-flanked300-clusters-map_positions.bed -fo ${controls}-5SS-flanked300-clusters-map_positions.fasta

# create a matrix of kmer and cluster coverage around each flanking region of regulated exons
python get_coverage-3SS.py ${exons}-3SS-flanked300-clusters-map_positions.fasta ${kmers} ${exons}-3SS-flanked300-clusters-map_positions
python get_coverage-5SS.py ${exons}-5SS-flanked300-clusters-map_positions.fasta ${kmers} ${exons}-5SS-flanked300-clusters-map_positions
python get_coverage-3SS.py ${controls}-3SS-flanked300-clusters-map_positions.fasta ${kmers} ${controls}-3SS-flanked300-clusters-map_positions
python get_coverage-5SS.py ${controls}-5SS-flanked300-clusters-map_positions.fasta ${kmers} ${controls}-5SS-flanked300-clusters-map_positions

# remove duplicates
sort ${exons}-3SS-flanked300-clusters-map_positions.csv | uniq > ${exons}-3SS-flanked300-clusters-map_positions-uniq.csv
sort ${exons}-5SS-flanked300-clusters-map_positions.csv | uniq > ${exons}-5SS-flanked300-clusters-map_positions-uniq.csv
sort ${controls}-3SS-flanked300-clusters-map_positions.csv | uniq > ${controls}-3SS-flanked300-clusters-map_positions-uniq.csv
sort ${controls}-5SS-flanked300-clusters-map_positions.csv | uniq > ${controls}-5SS-flanked300-clusters-map_positions-uniq.csv

# merge 3' and 5' splice site regions (in case that we have multiple clusters around one exon)
python merge_3SS_5SS-coveragee.py ${exons}-3SS-flanked300-clusters-map_positions-uniq.csv ${exons}-3SS-flanked300-clusters-map_positions-merged.csv
python merge_3SS_5SS-coveragee.py ${exons}-5SS-flanked300-clusters-map_positions-uniq.csv ${exons}-5SS-flanked300-clusters-map_positions-merged.csv
python merge_3SS_5SS-coveragee.py ${controls}-3SS-flanked300-clusters-map_positions-uniq.csv ${controls}-3SS-flanked300-clusters-map_positions-merged.csv
python merge_3SS_5SS-coveragee.py ${controls}-5SS-flanked300-clusters-map_positions-uniq.csv ${controls}-5SS-flanked300-clusters-map_positions-merged.csv

# draw a HeatMap
Rscript draw_HeatMap.R ${exons}-3SS-flanked300-clusters-map_positions-merged.csv ${exons}-5SS-flanked300-clusters-map_positions-merged.csv ${controls}-3SS-flanked300-clusters-map_positions-merged.csv ${controls}-5SS-flanked300-clusters-map_positions-merged.csv HeatMap-${clusters}-${exons}.pdf
