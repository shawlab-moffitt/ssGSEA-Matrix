


####----ssGSEA Score Generation----####

##--Input--##

Expression_File <- ""

GeneSet_File <- ""

Output_File <- ""



####----Packages----####

library(dplyr)
library(tidyr)
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

if (as.numeric(tools::file_path_sans_ext(packageVersion("GSVA"))) >= 1.5) {
  ssGSEA_param <- GSVA::ssgseaParam(expr,GeneSetList)
  ssGSEA <- GSVA::gsva(ssGSEA_param)
  ssGSEA <- as.data.frame(ssGSEA)
  ssGSEA$Genesets <- rownames(ssGSEA)
} else {
  ssGSEA <- gsva(expr,GeneSetList,method = "ssgsea", verbose = TRUE)
  ssGSEA <- as.data.frame(ssGSEA)
  ssGSEA$Genesets <- rownames(ssGSEA)
}
ssGSEA <- ssGSEA %>% relocate(Genesets)
write.table(as.data.frame(ssGSEA),Output_File, sep = '\t', row.names = F)
