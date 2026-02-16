# TRI Correlation

Folder contains dotmatrix plots generated via correlation of VX (day X) TRI
(titer response index) with gene counts or log2FoldChange gene counts.

Example filename I: `TRI_ELISA_V1_V1_HC_single_dotmatrix.png`  
depicts how V1 gene counts are correlated with the V1 TRI_ELISA.  

Example filename II: `TRI_ELISA_V4_V3_l2fc_count_corr_dotmatrix_combined_2.png`  
depicts how V1vsV3 log2FoldChanges in gene counts are correlated with the V4 TRI_ELISA.  

The difference between a `*_combined.png` and `*_combined_2.png` dotmatrix is the
method by which the pathways were chosen.  
For `*_combined.png` dotmatrices hierarchical clustering of pathways by gene
set similarity (Jaccard distances) was employed to reduce the number of pathways.  
For `*_combined_2.png` dotmatrices the top 10 highest and lowest NES (normalized enrichment score) pathways
were chosen from each diagnosis and merged.