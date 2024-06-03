**Sequencing samples:**  
Guide regions were amplified by low-cycle PCR from sample gDNA.  
Samples were sequenced by the UW-Madison Biotech Center Next Gen DNA Sequencing Core facility. PCR products were amplified with nested primers containing i5 and i7 indexes and Illumina TruSeq adapters, cleaned, quantified, pooled, and run on a Novaseq 6000 (150 bp paired-end reads). 

**Barcode counting info:**  
Guides were counted using a custom Python script designed by Ryan Ward to minimize noise and accurately count guides with overlapping targets. This script samples sequencing reads to identify barcode diversity, orientations, and offsets, and determines flanking sequences for correct barcode identification. Only counts for library barcodes were used in downstream analyses. This is far improved from bbtools for two reasons.  
1) "Fortuitous guides" can be identified and counted. Oligo synthesis is not perfect. Often, oligo libraries will contain random sequences that can, in fact, be active guides.  This script includes those.  
2) Guides that heavily overlap can be miscounted. Because this script samples to identify the appropriate flanking sequences and insertions, it ensures there aren't miscounts for guides which may be only 1-2 nucleotides offset from each other or for any contaminants/inappropriate insertions of guides.  
  
**Quality control**  
We observed instances of overselection, where too high of chemical concentrations selected for a underrepresented community of more fit strains that overtook the sample. This may be extremely useful in other experiments, but for the purposes of chemical genomics, this is too skewed to provide proper fitness information for library knockdowns. In a handful of cases, this overselection was highly replicable. To correct for this and any other noise, we removed samples with coefficients of determination (R^2) < 0.5 and if population diversity of nontargeting guides (Nb) was less than 10,000. In cases of overselection, the population diversity of nontargeting guides is extremely low, as most control guides in those samples were fully depleted due to antibiotic pressure.
  
**Relative fitness Scores**  
Log2 fold changes and confidence values were computed using EdgeR. Gene-level values were calculated by taking the median guide-level log2 fold-change for perfect match guides; confidence determined by computing Stouffer's p-value using guide-level FDR values (poolr R package).
