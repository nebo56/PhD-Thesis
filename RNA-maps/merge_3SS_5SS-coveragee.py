'''
Created on Aug 20, 2013

@author: Nejc

This script will merge columns from final .csv script where we have motif enrichment and clusters. In case
that we have multiple clusters in the same region we will merge them into a single row.
'''

import sys

# it will sum coverage of two arrays only if there is a difference
def sum_coverage(original_coverage, adding_coverage):
    for i in range(0, original_coverage.__len__()):
        if int(original_coverage[i]) < int(adding_coverage[i]):
            original_coverage[i] = int(adding_coverage[i])
    return original_coverage 

def merge(fin_fname, fout_fname):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    last_exon = None
    last_pos1 = None
    last_pos2 = None
    last_coverage = None
    line = fin.readline()
    while line:
        line = line.replace("-1.0","-1")
        col = line.rstrip('\n').rsplit(',')
        exon = col[0]
        pos1 = col[1]
        pos2 = col[2]
        str_coverage = col[3:]
        coverage = map(int, str_coverage)
        if last_exon != None:
            if exon != last_exon:
                fout.write(last_exon + ',' + last_pos1 + ',' + last_pos2 + ',' + str(last_coverage).replace(" ","").replace("[","").replace("]","") + '\n')
            else:
                if int(last_pos1) < int(pos1):
                    pos1 = last_pos1
                    pos2 = last_pos2
                coverage = sum_coverage(last_coverage, coverage)
        last_exon = exon
        last_pos1 = pos1
        last_pos2 = pos2
        last_coverage = coverage
        line = fin.readline()
    fout.write(exon + ',' + pos1 + ',' + pos2 + ',' + str(coverage).replace(" ","").replace("[","").replace("]","") + '\n')
   
    
if sys.argv.__len__() == 3:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    merge(fin_fname, fout_fname)
else:
    #print str(sys.argv.__len__())
    print "error:\t2 arguments are needed\n" + '\n' +"example:\t $ python merge_3SS_5SS-coveragee.py input_file.csv output_file.csv"
