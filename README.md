# Long-range-enh-gene-prediction-3layer
Probabilistic model and MCMC inference of long-range enhancer-gene regulation and associated TF modules
# Introduction
This repository contains code and predicted enhancer-gene links in each of 127 cell lines/tissues by integrative modelling combinatorial gene regulatory grammer. 
# Files and folders
`ENCODE_cell_index.csv`: metadata of cell lines, including ID, ENCODE sample name, epigenome mnemonic, standardized epigenome name and tissue.

`Code`: all scripts to predicted long-range enhancer-gene links. Including 8 steps of data preprocessing, 1 step for model, 4 steps for downstream validation and reformating.

   1. `1_generate_RNA_seq_profile.R` & `1_generate_DNase_profile.R`: generate gene expression and enhancer activity profile for each gene and enhancer.
   2. `2_generate_RNA_seq_matrix.R` & `2_generate_DNase_matrix.R`: reorganize results from step 1 into matrices for further use.
   3. `3_generate_motif_profile.R`: generate TF motif hits profile for each enhancer.
   4. `4_generate_potential_pair.R`: generate all potential enhancer-gene pool within certain distance window.
   5. `5_calculate_ep_corr.R`: calculate activity correlation for each potential enhancer-gene pair.
   6. `6_calculate_distance.R`: calculate distance between enhancers and genes for each pair.
   7. `7_extract_enh_TF_matrix.R`: reorganize results of step 3 into a matrix for further use.
   8. `8_prepare_input.R`: wrap up all the input data into a single .Rdata object.
   9. `9_integrative_model.R`: run integrative model.
   10. `10_generate_enh_promoter_frame.R`: get a data frame for all potential enhancer-gene links for later use.
   11. `11_ep_validation.R`: overlapping all potential enhancer-gene links with a certain gold-standard.
   12. `12_get_predicted_score.R`: extract predicted score of integrated model.
   13. `13_get_prediction.R`: get final predictions and write out as a table. 

`Prediction`: predicted enhancer-gene links for each of 127 cell lines/tissues. Imputed DNase and RNA-seq data in 127 cell lines/tissues are used as input. Predictions in this folder are selected based on the probability threshold at 0.85.

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

# Usage

A wrapper will be provided to run all scripts at once in the near feature. Currently the scripts can be run one by one in the order specified by the index.

## 1. `1_generate_RNA_seq_profile.R & 1_generate_DNase_profile.R`:
 
`Rscript 1_generate_RNA_seq_profile.R [path_to_exons/path_to_enhancer_coordinates] [index of cell types] [path to epigenomic datasets]`
### Notes:
`path_to_exons/path_to_enhancer_coordinates`: a cohort of exons/consensus enhancers you would like to use for enhancer-gene link prediction. We suggest to use exons from GENCODE V19 and enhancers from Roadmap. See `Input data description` below for downloading URLs.
`index of cell types`: from 1 to N (number of cell types of the input profile)
`path_to epigenimic datasets`: directory where you store the DNase-seq and RNA-seqs in .bedGraph format.

## 2. `2_generate_RNA_seq_matrix.R`: 
`Rscript 2_generate_RNA_seq_matrix.R`
###Notes:
By default the gene annotation is GENCODE V19, which is stored under `gene_annotation`. You can switch to a different gene model by overwritting the original file with the same file name and format. The format is descriped in the `Input data description`.

## 3. `3_generate_motif_profile.R`:
`Rscript 3_generate_motif_profile.R [path_to_enhancer_coordinates] [chromosome_index]`
###Notes:
The script will overlap motif with enhancers within a specific chromosome, will need to run the scripts for all 23 chromosomes seperately.

## 4. `4_generate_potential_pair.R`: 
`Rscript 4_generate_potential_pair.R [path_to_enhancer_coordinates] [path_to_promoter_locations]`
###Notes:
The default promoter location is defined by extending TSS of each protein coding gene by +/- 1kb. The TSS is derived from GENCODE V19 annotation. Different gene model can be used by overwriting this file.

## 5. `5_calculate_ep_corr.R`: 
`Rscripts 5_calculate_ep_corr.R`

## 6. `6_calculate_distance.R`:
`Rscripts 6_calculate_distance.R`

## 7. `7_extract_enh_TF_matrix.R`:
`Rscript 7_extract_enh_TF_matrix.R [chromosome_index]`

## 8. `8_prepare_input.R`:
`Rscript 8_prepare_input.R [index_of_cell_types] [max_distance_between_enhancers_and_promoters]`
###Notes:
The maximum distance should between 0 and 2.000,000 (2MB).

##9. `9_integrative_model.R`:
`Rscript 9_integrative_model.R [Initial_weight_of_TF_profile_for_ep_links] [Initial_weight_of_TF_profile_for_gene_group] [C] [number_of_gene_groups] [output_path] [min_tolerance_of_KL_divergence] [input_data_path] [d_m] [Randomize_methods]`
###Notes:
To be updated.

## 10.`10_generate_enh_promoter_frame.R`:
`Rscript 10_generate_enh_promoter_frame.R`

## 11. `11_ep_validation.R`:
`Rscript 11_ep_validation.R [index_of_block: 1 to 500] [path_to_gold_standards]`
###Notes:
Due to computational issue, all potential enhancer-gene links are divided into 500 blocks. The script can do overlapping for the block specified by the index.

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

