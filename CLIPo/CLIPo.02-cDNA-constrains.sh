#!/bin/bash -l

cDNAs=$1
path=/SAN/neuroscience/TE/Nejc/CLIPo-all/

gunzip ${path}${cDNAs}.gz
python ${path}scripts/separate_cDNAs_40nt_length.py ${path}${cDNAs} tmp

complete=tmp.40less.bed
incomplete=tmp.40more.bed

cat ${complete} | grep '-' > tmp.minus.bed
# separate plus and minus strand
python ${path}scripts/separate_cDNAs_by_strand.py ${complete} tmp

# sort by cDNA end positions
sort -k3,3n tmp.plus.bed > tmp.plus.cDNAend.sort.bed
sort -k2,2n tmp.minus.bed > tmp.minus.cDNAend.sort.bed
cat tmp.plus.cDNAend.sort.bed tmp.minus.cDNAend.sort.bed > tmp.cDNAend.sort.bed

# get score for cDNA ends (how many cDNAs ends are on the same position)
python ${path}scripts/get_cDNAend_score-count.py tmp.cDNAend.sort.bed tmp.cDNAend.sort-score.bed

# add incomplete cDNAs to our dataset !!! for eCLIP change $5 to $4 and the other way around
cat ${incomplete} | awk '{print $1 "\t" $2 "\t" $3 "\t""\t""\t" -1 "\t" $5}' > tmp.incomplete.bed
cat tmp.cDNAend.sort-score.bed tmp.incomplete.bed > tmp.cDNAend.sort-score-incomplete.bed

#separate by strands
python ${path}scripts/separate_cDNAs_by_strand.py tmp.cDNAend.sort-score-incomplete.bed tmp.cDNAend.sort-score

#sort by cDNA start positions
sort -k2,2n tmp.cDNAend.sort-score.plus.bed > tmp.cDNAstart.sort.plus.bed
sort -k3,3n tmp.cDNAend.sort-score.minus.bed > tmp.cDNAstart.sort.minus.bed
cat tmp.cDNAstart.sort.plus.bed tmp.cDNAstart.sort.minus.bed > tmp.cDNAstart.sort.bed

# get score for cDNA starts (how many cDNAs starts are on the same position)
python ${path}scripts/get_cDNAstart_score-count.py tmp.cDNAstart.sort.bed ${path}${cDNAs}-cDNAscore.bed

# separate by strand and sort plus strand by starts and minus strand by ends
#python separate_cDNAs_by_strand.py ${cDNAs}-cDNAscore.bed tmp.cDNAscore
#sort -k1,1 -k2,2n -k6,6 tmp.cDNAscore.plus.bed > tmp.cDNAscore.plus.sorted.bed
#sort -k1,1 -k3,3n -k6,6 tmp.cDNAscore.minus.bed > tmp.cDNAscore.minus.sorted.bed
#cat tmp.cDNAscore.plus.sorted.bed tmp.cDNAscore.minus.sorted.bed > tmp.cDNAscore.both.sorted.bed

# separate if cDNAs start socre is >= cDNAs end score
python ${path}scripts/separate-cDNAs-score-report.py ${path}${cDNAs}-cDNAscore.bed > ${path}${cDNAs}-02.CLIPo-cDNA-constrains-REPORT.txt

rm tmp.*
