


####----ssGSEA Score Generation----####

##--Input--##

Expression_File <- "~/R/data/GrueberNeault_data/Gruber_Pediatric_hememalignancies/Human_Clinical_AML_BALL_HM.combined.fpkm.geneNameOnly.max.reorder.txt"

GeneSet_File <- "~/R/data/GrueberNeault_data/Gruber_Pediatric_hememalignancies/TRRUST_Transcription_Factors_2019.gmt"

Output_File <- "~/R/data/GrueberNeault_data/Gruber_Pediatric_hememalignancies/Gruber_Pediatric_hememalignancies_TRRUST_Transcrpt_Factors_ssGSEA.txt"



####----Packages----####

library(dplyr)
library(readr)
library(GSVA)
library(clusterProfiler)



##--Read In Files--##

## Expression
expr <- as.data.frame(read_delim(Expression_File,delim = '\t', col_names = T))
expr <- expr %>%
  drop_na() %>%
  as.data.frame()
# Check that expression data is numeric
isChar <- unname(which(sapply(expr, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr[isChar] <- sapply(expr[isChar],as.numeric)
}
colnames(expr)[1] <- "Gene"
# Remove Duplicate genes
if (TRUE %in% duplicated(expr[,1])) {
  expr <- expr %>%
    group_by(Gene) %>%
    summarise_all(max) %>%
    as.data.frame()
}
row.names(expr) <- expr[,1]
expr <- as.matrix(expr[,-1])

## Gene Set
if (tools::file_ext(GeneSet_File) == "gmt") {
  GeneSet <- read.gmt(GeneSet_File)
  colnames(GeneSet) <- c("term","gene")
  GeneSetList <- list()
  for (j in unique(GeneSet[,1])){
    GeneSetList[[j]] <- GeneSet[GeneSet[,1] == j,]$gene
  }
}
# If user loads TSV/TXT file
if (tools::file_ext(GeneSet_File) == "tsv" || tools::file_ext(GeneSet_File) == "txt") {
  GeneSet <- read.delim(GeneSet_File, header = T, sep = '\t')
  colnames(GeneSet) <- c("term","gene")
  GeneSetList <- list()
  for (j in unique(GeneSet[,1])){
    GeneSetList[[j]] <- GeneSet[GeneSet[,1] == j,]$gene
  }
}
# If user loads RData list file
if (tools::file_ext(GeneSet_File) == "RData") {
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  GeneSetList <- loadRData(GeneSet_File)
}



ssgsea <- gsva(expr,GeneSetList,method = "ssgsea")   #Outputs sample names as column names
ssgsea <- as.data.frame(t(ssgsea))                   #Transpose to make column names gs names
ssgsea$SampleName <- rownames(ssgsea)                #Change gs names to row names
ssgsea <- ssgsea %>%                                 #Make rownames first column to output to file
  relocate(SampleName)
write_delim(as.data.frame(ssgsea),Output_File, delim = '\t')

