# IgG Correlation

Folder contains dotmatrix plots generated via correlation of VX (day X) IgG titers
(as measured via ECL-ELISA) with gene counts or log2FoldChange gene counts.

Example filename I: `IgG_V1_V1_HC_single_dotmatrix.png`  
depicts how V1 gene counts are correlated with V1 IgG levels.  

Example filename II: `IgG_V4_V3_l2fc_count_corr_dotmatrix_combined_2.png`  
depicts how V1vsV3 log2FoldChanges in gene counts are correlated with V4 IgG levels.  

The difference between a `*_combined.png` and `*_combined_2.png` dotmatrix is the
method by which the pathways were chosen.  
For `*_combined.png` dotmatrices hierarchical clustering of pathways by gene
set similarity (Jaccard distances) was employed to reduce the number of pathways.  
For `*_combined_2.png` dotmatrices the top 10 highest and lowest NES (normalized enrichment score) pathways
were chosen from each diagnosis and merged.
