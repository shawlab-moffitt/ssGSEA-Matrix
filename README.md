# ssGSEA-Matrix


# Scripts

## ssGSEA_Scoring.R
This is an R script that will perform a ssGSEA calculation on a given expression matrix with a given geneset file.
### Requirments
|  |  |  |
| --- | --- | --- |
| readr_2.1.2 | dplyr_1.0.9 | tidyr_1.1.3 |
| GSVA_1.40.1 | clusterProfiler_4.0.5 |  |
### Input
* **Expression_File**: Expression matrix with the first column consisting of feature names and the subsiquent columns with a header of sample names and the values contain numeric feature scores.
* **GeneSet_File**: A gene set file in gmt, two-column tsv/txt, or an RData list format.
  * **gmt**: First column: geneset names, second column: description, subsiquent columns genes within the geneset. Example found [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
  * **two-column tsv/txt**: The first column being the gene set name repeating for every gene symbol that would be placed in the second column. Example found [here](https://raw.githubusercontent.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/main/GeneSets/CellMarker_gsNsym_HS.tsv)
  * **RData list**:A named list of gene sets and genes. Example found [here](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/tree/main/3-PATH-SURVEYOR-Pipeline/GeneSets/CellMarker). A script to generate this list is provided here: [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/PATH-SURVEYOR-Suite/blob/main/3-PATH-SURVEYOR-Pipeline/GeneSets/GeneSetRDataListGen.R). 
* **Output_File**: Desired name of output file
### Output
Tab-delimited matrix file where the first column consists of geneset names and the subsiquent columns are each sample with the numeric value of their ssGSEA score per gene set.
