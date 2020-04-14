# Long-range-enh-gene-prediction-3layer
Probabilistic model and MCMC inference of long-range enhancer-gene regulation and associated TF modules
# Introduction
This repository contains predicted enhancer-gene links in each of 127 cell lines/tissues by integrative modelling combinatorial gene regulatory grammer. 
# Files and folders
`ENCODE_cell_index.csv`: metadata of cell lines, including ID, ENCODE sample name, epigenome mnemonic, standardized epigenome name and tissue.

`Prediction`: predicted enhancer-gene links for each of 127 cell lines/tissues. 

1. Two versions are provided in corresponding sub-directory: 
   1. 0.85: predictions in this folder are selected based on the probability threshold at 0.85. The threshold corresponds to  p-value=0.05 in the permutation test. Roughly there will be ~100,000 enhancer-gene links in each file.
   1. distance_0.8_ChIA_PET: predictions in this folder are selected based on both distance and probability. Specifically, top-ranked predictions following **the same distance distribution as the experimental ChIA-PET data in K562 (GSE59395)** are selected. All selected predictions had a minimum probability of 0.8. Roughly there are ~50,000 enhancer-gene links in each file.
2. `Cell_index.txt`: prediction for each of 127 cell lines. The detailed information of cell lines/tissues can be retrieved using ENCODE_cell_index.csv. The file is organized as below:

Column index | Column description | Column type
------------ | ------------- | -------------
1 | Enhancer chromosome name | character
2 | Enhancer start | int
3 | enhancer end | int
4 | Target promoter name | character
5 | Target promoter chromosome name | character
6 | Target promoter start | int
7 | Target promoter end | int
8 | distance | int
9 | correlation | float
10 | Probability of linking | float

# Input data description:

The integrative model will predict enhancer-gene links based on: (1) activity correlation between enhancers and genes across 127 cell lines/tissues; (2) distance between enhancers and genes; (3) TF motif hits within enhancer regions. Detailed description of all the input data are provided below:

1. *Define consensus enhancer coordinates across 127 cell lines/tissues*

   The consensus enhacer coordinates were downloaded from the website https://personal.broadinstitute.org/meuleman/reg2map/HoneyBadger2-intersect_release/, with the version of DNaseI regions selected with -log10(p) >= 10. These regions are delineated using observed DNaseI data across 53 epigenomes and annotated with both the 5-mark 15 states model based on observed data as well as the 12 mark 25-state model based on imputed data, across 127 epigenomes. Totally it marks 474,004 putative consensus enhancer regions across 127 cell types.

2. *Generate enhancer activity profile across 127 cell lines/tissues.*

   The enhancer activity of each enhancer is defined as the average DNase signal within the enhancer region in each of 127 cell lines/tissues. The imputed whole genome-wide DNase signal tracks were provided by the Roadmap project (https://egg2.wustl.edu/roadmap/web_portal/imputed.html#imp_sig). The data was downloaded in ‘.bigwig’ format and converted to ‘.bedGraph’ format with bedTools. To calculate the DNase signal for each enhancer, genomic bins from ‘.bedGraph’ files were overlapped with enhancers. The averaged signal across overlapping bins was used as the enhancer activity. The enhancer activity profile was organized into a matrix, where each row represents an enhancer and each column represents one cell types.

3. *Generate gene expression profile across 127 cell lines/tissues.*

    The Gencode V19 was used for gene annotation. Only protein-coding genes on chr1 to chrX were used for subsequent analyses (~20,300 genes).  The promoter was defined as the +/-1kb region centered at the TSS of genes. 
    The imputed RNA-seq data for all 127 cell lines/tissues were provided by the Roadmap (https://egg2.wustl.edu/roadmap/web_portal/imputed.html#imp_sig). The downloaded ‘.bigwig’ format was converted to ‘.bedGraph’ format with bedTools. To quantify the expression level of each gene, the average RNA-seq signal (log(RPKM)) within all exons for each gene was calculated. Exons are extracted from the Gencode V19 annotation file. Since the ENSEMBL annotation only marks exons for ~7000 out of ~20000 protein coding genes, both HAVANA and ENSEMBL annotations from the original ‘gtf’ file were used. To remove overlapped and duplicated exons coming from different resources (HAVANA or ENSEMBL), exons for each gene were merged and the resulting non-overlapping exons were used for later calculation. The calculation of averaged signals of all exons for each gene was the same as that stated in section 2. 
The results of both methods were organized into a matrix, where each row is one gene and each column is one cell type.

4. *Generate potential enhancer-promoter links and distances*

    For each promoter, links between the TSS and all enhancers that are within +/- 2MB distance window centered at the TSS were considered as potential enhancer-promoter links. The distance between an enhancer and the TSS was defined as the distance between the center of the enhancer and the TSS.

5. *Generate correlations between enhancers and promoters across 127 cell lines/tissues*

    The DNase profile of enhancers that were generated in section 2 and the RNA-seq profile that were generated in 3 were quantile normalized across 127 cell lines/tissues (by column). For each enhancer-promoter link, the pearson correlation was calculated based on the normalized enhancer’s DNase vector and promoter’s RNA-seq signal vector.

6. *Generate the TF motif hit profile for each enhancer*

    The genome-wide motif annotation data was downloaded from http://compbio.mit.edu/encode-motifs/matches-with-controls-0.3.txt.gz. Controls (indicated by “_C”) were removed. This resulted in ~13.6M motif hits for ~500 TFs. These motif hits were overlapped with enhancers. In each of 127 cell lines/tissues, only expressed TFs (log(RPKM)>1) were retained and used as input of the integrative model. The results were organized into a matrix, where each row is one enhancer, each column is one TF and each entry shows the number of motifs of each TFs in each enhancer. For each of the 127 cell lines, only columns corresponding to expressed TFs (RNA-seq signal: imputed log(RPKM) >1, based on exon- or promoter-metrics) were kept and used as input of the model.

