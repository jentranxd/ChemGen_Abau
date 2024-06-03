Once again, much of this code takes inputs direct from variables created in the EdgeR script

For pathway comparisons, more information can be found here: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html  
  
For structural comparisons, more information can be found here: https://chemminetools.ucr.edu/help/ 

**Brief description:**  
After edgeR was used to import, organize, filter, and normalize counts, the limma package with voom method was used to perform gene set testing using functional groups predicted by STRING-db (custom proteome; Organism ID: STRG0060QIE). Chemical effects on pathways were scored for direction and significance, then transformed into a correlation matrix. Chemicals were clustered using hierarchical clustering with the Ward method, and a dendrogram visualized the clustering based on pathway impacts.
SMILES strings for each chemical were parsed to generate molecular fingerprints. For combination antibiotics, a representative structure was selected. Tanimoto coefficients were calculated to measure chemical similarity and generate a similarity matrix, which was then clustered using the Ward method to produce a chemical structure dendrogram.
To assess the relationships between chemical pairs, each dendrogram was iteratively cut at all k-values (1-45), and pairs were scored if they clustered together. The final association scores are the sums across all k-values. Differences in association scores between pathway-effect and structure were plotted to highlight disparities. Additionally, the sums of association scores for pathway-effect and structure were plotted to highlight similarities.

