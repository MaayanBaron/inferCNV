# inferCNV

## This script infers Copy Number variation (CNV) from scRNA-seq data inspired by Tirosh et al. 2014 (just in Matlab). The output is a figure showing the CNV profile for each cells (every row is a cell) 
 
Dependencies: MAGIC as an imputution method https://github.com/KrishnaswamyLab/MAGIC

## INPUT ARGUMENTS
   ---------------
   A_fil - a filtered expression matrix (cells as columns, genes as rows)
   clustLabel - a numeric vector indicating the cluster for each cell in A
   (cancer, immune, stroma...)
   clusNames - a string vector indicating the cluster names for each of the
   values of clustLabels
   gene_names - all of the gene names of A_fil
   gene_loc - a matrix in the length of gene_names containg 3 columns:
   start location, end location and chr number
   gene_name_ord - all of the gene names of gene loc
   sample_ID_fil - sample id (for multiple batches)
   set of genes to remove - if wanted, you can remove specific genes from
   the analysis (for example: differentially expressed genes)

## OUTPUT ARGUMENTS
   ----------------
   The results of this script is infered CNV from the scRNA-seq plotted 
   on a graph

