# PhD thesis
# Chapter 3: Assessing potential biases in protein-RNA binding site assignment with iCLIP

# Trimming of adapter sequences
Before mapping the cDNAs, we removed random barcodes and trimmed the 3' Solexa adapter sequence. Adapter sequences were trimmed by FASTX-Toolkit 0.0.13 adapter removal software, using the following parameters: -Q 33 -a AGATCGGAAG -c -n -l 26. For reads that did not contain parts of the adapter sequence, "-C" parameter was used, and these were analyzed separately.

# Genome mapping of PTBP1 and U2AF65 iCLIP
We used UCSC hg19/GRCh37 genome assembly and bowtie2 2.1 alignment software with default settings accepting uniquely mapped cDNAs to a single genomic position and allowing maximum of 2 missmatches. After mapping, cDNAs with the same random barcode that mapped to the same starting position on the genome were considered to result from PCR amplification and were collapsed to a single cDNA.

# Transcriptome mapping of eIF4A3 iCLIP
For mapping, we compiled a set of representative mRNA sequences from BioMart Ensembl Genes 79, where we used the longest mRNA sequence available for each gene. We mapped to these mRNAs with Bowtie2.1 alignment software, allowing 2 mismatches. After mapping, cDNAs with the same random barcode that mapped to the same starting position on an mRNA were considered to result from PCR amplification and were collapsed to a single cDNA.

# Classification of cDNA length
Only cDNAs that mapped to a unique genomic position were evaluated. These were separated into cDNAs that were <30 nt, 30-34 nt, 35-39 nt or >40 nt long after trimming. 

# Definition of crosslink-associated motifs
We reasoned that sequence motifs enriched directly at the starts of the control eCLIP cDNAs might uncover preferences of UV crosslinking, since they are thought to represent a mixture of crosslink sites for many different RBPs, and thus they should not reflect sequence specificity of any specific RBP. We therefore examined occurrence of tetramers that overlapped with the nucleotide preceding the cDNA-starts in PTBP1 control iCLIP (position -1) in comparison with the ones overlapping with the 10th nucleotide preceding the cDNA-starts (position -10). Tetramers that are enriched over 1.5 fold at position -1 compared to -10 include. We excluded the TTTT tetramer from further analyses, since it is often part of longer tracts of Ts, and therefore its inclusion decreases the resolution of analysis. Thus, TTTG, TTTC, TTGG, TTTA, ATTG, ATTT, TCGT, TTGA, TTCT and CTTT were used for all analyses of crosslink-associated motifs.

# Definition of PTBP-target motifs
To identify the motifs bound by PTBP1, we searched for pentamers enriched in the region [-10..10] around the cDNA-start peaks identified in each crosslink cluster defined by PTBP1-iCLIP2.  69 pentamers had enrichment z-score > 299 and were used as PTBP1-target pentamers for further analyses. Their sequences are: TCTTT, CTTTC, TCTTC, CTTCT, TCTCT, CTCTC, TTTCT, TTCTC, TTCTT, TTTTC, TCCTT, CTCTT, ATTTC, TTCCT, CTTCC, TTTCC, CCTTT, CTTTT, CCTTC, TCTGT, TTCTG, TCCTC, CTTCA, ATCTT, TGTCT, TCTGC, CTCCT, CCTCT, GTCTT, TCTAT, TCTCC, ATTCC, TTCTA, CTTTG, TATCT, ACTTC, TTATC, CTTAT, CTATT, TTCAT, TTCCA, TCTTG, TTGTC, TTGCT, CTCTA, CTCTG, TATTT, TCCCT, TCATT, TTCCC, CATTT, ATTCT, TTTAC, GTTCT, CTATC, TCATC, CTTTA, TGTTC, TATTC, CATCT, TACTT, CTGTT, CTTGC, ACCTT, TTTCA, TTTGT, TGTTT, CTTGT, ACTTT. All of these pentamers are enriched in pyrimidines, in agreement with the known preference of PTBP1 for UC-rich binding motifs (see "Y-rich-pentamers.tab" in the Chapter4/RNA-maps folder).

# Normalisation of data for drawing of density graphs
All normalisations were performed in R (version 3.1.0) together with “ggplot2” and “smoother” package for final graphical output.

- For analysis of EIFA3 iCLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of exon-exon junctions in mRNAs. All exon-exon junctions within coding regions were taken into account, apart from those that junctioned to the first or last exon in the mRNA. The number of cDNAs starting or ending at each position on the graph was normalised by the number of all cDNAs mapped to representative mRNAs, the mRNA length and the number of examined exon-exon junction positions in the following way: 
RNAmap[n] = ((cDNAs[n] / sum(cDNAs)) * length(mRNAs) / count(exon_exon_junctions)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 3-nucleotide smoothing window. 

- For analysis of PTBP1, U2AF65 iCLIP and CLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of its binding sites. We defined the binding sites in two ways: either by using the position of Y-tracts or extended crosslink clusters in introns. All introns in protein coding genes were taken into account. Due to highly variable abundance of intronic RNAs (and occasional presence of highly abundant non-coding transcripts, such as snoRNAs), we first divided counts at each binding site by dividing them by the count at the MaxCount. To define MaxCount, we examined the region of the binding site, as well as 120 nt 5' and 3' of the binding site, to find the nucleotide with the largest count of cDNA count starts or cDNA ends (according to whether starts or ends were plotted on the graph). This avoids from the highly abundant intronic snoRNAs or other abundant introns from dominating the results. The MaxCount-normalised counts of cDNAs starting or ending at each position on the graph was then normalised by the density of all MaxCount-normalised cDNA count/nt starting in the region 50-100 nt downstream of the binding site in the following way: 
RNAmap[n] = (MaxCount-normalised cDNAs[n]) / MaxCount-normalised cDNA density(outside the binding site)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 10-nucleotide smoothing window. 

# Assignment of the cDNA-end and cDNA-start peaks in eIF4A3 iCLIP
For cDNA-end peak assignment in eIFA3 iCLIP data, we used exons longer then 100 nt that were in the top 50% of the distribution of exons based on cDNA coverage. This ensured that sufficient cDNAs were available for assignment of the putative binding sites. We then summarised all cDNA-end positions in the region -20 to +25 around exon-exon junctions and selected the position with the maximum cDNA count as the ‘cDNA-end peak’.

# Analysis of pairing probability
Computational prediction of the secondary structure around the cDNA-end peaks was performed using the RNAfold program with the default parameters.

# Analysis of cDNA transitions
Density of U>T transitions across cDNAs was performed by using the samtools software with the the following parameters: samtools calmd –u –u genomic.fasta input_BAM > BAM-with_transitions. This pipeline replaces BAM format mapped cDNA sequences with transitions relative to genomic reference. This step was performed with a custom python script (available on github repository) that returns a density array of U>T transitions for cDNAs that are shorter then 40 nts. For the final visualisation of density graphs we used the same approach as for all other density figures. 

# Identification of crosslink clusters
The crosslink clusters were identified by False Discovery Rate peak finding algorithm from iCount (https://github.com/tomazc/iCount), by considering all crosslink sites that were significant with a FDR<0.05 at a maximum half windows spacing of 3 nt between crosslink sites. Than, the significant peaks were merged into final clusters by distance of 3 nt.


# Chapter 4: CLIPo: a tool to identify the features underlying protein-RNA interactions from CLIP data

# Visualisation of RNA-maps
This method takes into account predefined crosslink cluster positions, exonic positions that are regulated by the RBP of our interest, together with unregulated control exons and motif sequences. The pipeline can also be used with non-regulated control exons or with out motifs (RBP kmers), which are used for visualisation purposes. In the first step, employs bedtools intersect function to select all the neighbouring clusters in 300 nt flanking region around regulated exons, and then sort them by the distance from exon starts. In the second step, it extracts the genomic sequences around the splice site regions and creates a matrix based on cluster positions and motif enrichment (Y-rich motifs were used for PTBP1 and hnRNPC). This part of the process is written in Python and it takes each nucleotide position around regulated exons to set a value based on their cluster position and motif enrichment: -1 value is for exon start or end position, value 1 is each position that overlaps with any of the motifs, 2 is for cluster positions, 3 is for the motif coverage that is inside of a cluster region and every other position is set to 0 value. These values are stored as a matrix in CSV format and are then visualised with R script, where the matrix is plotted as a heat map, a density plot of cluster enrichment and a table. In the heatmap, every row represents a regulated exon with the cluster position and the motif coverage in surrounding region 300 nt upstream and downstream from the regulated exons and 50 nt inside of an exon. The second plot is a density plot of clusters enrichment compared to the enrichment in control exons. The last result of this pipeline is a table with crosslink cluster enrichments with distances and ratios between control exons and regulated exons in 3’ splice site region, inside of exons and 5’ splice site region.
Analysed RBPs:

# CLIPo analysis
Results from the CLIPo table in Chapter 4 were calculated with the following analysis.

- Data complexity: 
Library size of uniquely mapped cDNAs after PCR duplicate removal.
- cDNA constrains: 
I only focus on cDNAs that had a full length sequenced by looking at those which are less than 40 nt long.

- Length constraints: 
I used a sliding window approach to detect the most enriched cDNA length density in a 10 nt window frame. This value tells us if there are some strong cDNA length constrains such as narrow cDNA lengths group across the library.
- Sequence constraints at cDNA-ends: 
Ratio of top 10 tetramers that are positioned around cDNA-ends (from less than 40 nts long cDNAs) compared to the top 10 tetramers from central region between -15 to -5 nt upstream from cDNA end.
R=(∑top10 at cDNA.end)∗6(∑top10 at cDNA.centre)(3)
(3)R=(∑top10 at cDNA.end)∗6(∑top10 at cDNA.centre)

- Structure constraints at cDNA-ends: 
Ratio of single strandedness around cDNA-ends (from less than 40 nt long cDNAs) compared to the central region between -15 to -5 nt upstream from cDNA end. Single strandedness was measured same as above in a 30 nt surrounding region.
- Number of crosslink clusters: 
The number of crosslink clusters was obtained with iCount peak calling tool. The clusters were identified as described earlier, by considering all crosslink sites that were significant with a FDR < 0.05 but with a maximum spacing of 15 nt between crosslink sites.
- Percentage of cDNA-starts in the clusters: 
Measured percentage of cDNA-starts that are inside of identified crosslink clusters.
- Motif enrichment inside the clusters: 
Enrichment of tetramer coverage between identified clusters and 300 nt control regions downstream from the clusters. Top 10 tetramers were selected from those identified around cDNA-starts in 10 nt surrounding region.

# Chapter 5: Assignment of RNA binding sites for higher-order proteins complexes

# Identification of cDNA-start peaks and tetramer enrichment
I processed each mapped PTBP1 iCLIP, eCLIP and mock-eCLIP dataset with the iCount pipeline to define crosslink clusters with 3nt spanning window, 20 nt cluster merging and 0.05 FDR threshold. I only selected cDNAs that were inside of those clusters and then I selected position with the highest cDNA-start count for each cluster and defined it as a cDNA-start peak for further analysis. Next, I ranked all tetramers that were enriched in 20 nt flanking region around the maximum peaks. The enrichment of each tetramer was measured in comparison with the control frequency of tetramers from non-overlapping region of 200 to 300 nt downstream from cDNA-start peaks. I used the same peaks with the same surrounding region and controls to measure the enrichment of pairing probability using RNAfold software and a python script as described before. For the correlation between tetramer enrichments I used Pearson correlation. I calculated the individual upper and lower quartile of cDNA-start peaks for the most common tetramers and used them for further analysis. The same conditions were used for the pairing probability analysis.

# All scripts
Chapter 3: 
 - mapping_to_genome-PTBP1_and_U2AF65_iCLIP-pipeline.sh (PTBP1 and U2AF65 mapping pipeline to genome)
 - mapping_to_genome-PTBP1_and_U2AF65_CLIP-pipeline.sh (PTBP1 and U2AF65 mapping pipeline to genome)
 - mapping_to_transcriptome-eIFA3_iCLIP-pipeline.sh (eIFA3 mapping pipeline to transcriptome)
 - mapping_to_transcriptome-eIFA3_CLIP-pipeline.sh (eIFA3 mapping pipeline to transcriptome)
 - normalisation_and_density_graph-eIFA4.R (normalisation and drowing of eIFA3)
 - normalisation_and_density_graph-PTBP1-U2AF65.sh (normalisation and drowing of PTBP1 and U2AF65)
 - main_get_cDNAstart-end_peaks.sh (eIFA3 selecting top 1000 transcripts and reporiting/drawing around cDNA start and cDNA end peaks)
 - get_cross-link_clusters.sh (get cross-link clusters from mapped cDNAs with iCount - peak calling function)
 - make_HeatMap.sh (PTBP1 motif heatmaps for grouped clusters)

Chapter 4: CLIPo
- CLIPo.01.mapping_to_genome-CLIP.sh (mapping pipeline for CLIP data)
- CLIPo.01.mapping_to_genome-iCLIP.sh (mapping pipeline for iCLIP data)
- CLIPo.01.mapping_to_genome-eCLIP.sh (mapping pipeline for eCLIP data)
- CLIPo.01.mapping_to_genome-irCLIP.sh (mapping pipeline for irCLIP data)
- CLIPo.02-cDNA-end_constraints.sh (report fo cDNA-end constraints)
- CLIPo.02-cDNA-end_structure_constraints.sh (report fo cDNA-end structure constraints)
- CLIPo.03-kmer-constraints.sh (report fo cDNA-end/end kmer contraints)
- CLIPo.04-identify_crosslink_clusters-iCount.sh (pipeline for crosslink cluster identification with iCount peak calling tool)
- CLIPo.05-clusters-specificity.sh (kmer and cDNA-start coverage across all crosslink clusters)

Chapter 4: RNA-maps
- get3SS.py (get exons 3' splice site)
- get5SS.py (get exons 5' splice site)
- flankBEDpositionsCustom.py (flank splice site positions)
- get_coverage-3SS.py (get coverage of kmers and clusters around exons 3' splice sites)
- get_coverage-5SS.py (get coverage of kmers and clusters around exons 5' splice sites)
- get_flanked_cluster-positions.py (get positions of clusters from flanked splice sites)
- merge_3SS_5SS-coverage.py (merge 3' and 5' flanking splice sites)
- draw_HeatMap.R (draw a final heatmap)
- make_RNA-maps.sh (main script)

Chapter 5: Branch Points
- BEDsum.py (add count based on the number of overlapping cDNA-start positions)
- BranchPoint-detection.sh (main script)
- flankBEDpositions.py (flank splice site positions)
- get_branch_point_candidates.py (filter cDNA sequences that ends with AG)
- set_branch_point_position.py (set 1 nt position upstream from cDNA-start as branch point position)
- swap_barcodes_to_header.py (move random barcode UMI from the iCLIP fasta sequence to its header)
- transition_ratio.py (count a genomic transitions on first nucleotide from each read. Input file must be in SAM format)
- trimSAM.py (trim cDNAs that contain a genomic A mutation at the first nucleotide in the SAM sequence file)


# Common scripts
 - other scripts (kmer finder, flanking BED positions, density of deletions across all cDNAs, density of C to T transitions)
