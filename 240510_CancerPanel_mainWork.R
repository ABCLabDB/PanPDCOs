####################################
## Analysis for cancer panel-seq  ##
## 2024-05-10                     ##
## Minseok Seo                    ##
####################################

## Load packages
setRepositories(ind = 1:7)

library(data.table)
library(edgeR)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
library(stringr)
library(corrplot)
library(pheatmap)
library(readxl)
library(tidyverse)

## Set Working Dir.
WORK_DIR <- ""
DATA_DIR <- "./0.Data/2.CancerPanelSeq"

setwd(DATA_DIR)

## Load raw variant data derived from Macrogen
FILENAME <- "AnnotatedVariants.xlsx"
NUMSAMPLES <- 11

newLabel <- c("128_FNB_PDCO_S",
              "128_ERCP_PDCO_S",
              "128_OP_FFPE",
              "124_ERCP_PDCO_S",
              "124_FNB_PDCO_S",
              "124_OP_FFPE",
              "134_FNB_PDCO_S",
              "134_ERCP_PDCO_S",
              "134_OP_FFPE",
              "116_OP_PDCO_S",
              "116_OP_FFPE")

data <- list()

for(i in 1:NUMSAMPLES){
  data[[i]] <- read_excel(FILENAME, sheet = i)
}

names(data) <- excel_sheets(FILENAME)

data.frame(names(data))

## Extract required info.
for(i in 1:NUMSAMPLES){
  data[[i]] <- data.frame(ID = paste0(data[[i]]$CHROM, ":", data[[i]]$POS),
                          CHROM = str_remove(data[[i]]$CHROM, "chr"),
                          POS = as.numeric(data[[i]]$POS),
                          FILTER = data[[i]]$FILTER,
                          AnnoGene = data[[i]]$Gene)
}

## Remove duplicated variants
for(i in 1:NUMSAMPLES){
  data[[i]] <- data[[i]][!duplicated(data[[i]]$ID),]
}

## Filtering artifacts
selectFilter <- c("germline;panel_of_normals", "PASS", "haplotype;panel_of_normals", "haplotype", "germline;haplotype;panel_of_normals", "germline", "germline;haplotype", "panel_of_normals")

for(i in 1:NUMSAMPLES){
  data[[i]] <- data[[i]][which(data[[i]]$FILTER %in% selectFilter),]
  print(dim(data[[i]]))
}


## Remove non somatic mutations (Option)
#for(i in 1:NUMSAMPLES){
#  data[[i]] <- data[[i]][which(data[[1]]$FILTER == "PASS"),]
#}

## Full outer join
data <- data %>% 
  reduce(full_join, by = c('ID', 'CHROM', 'POS', 'FILTER'))

data$CHROM[which(data$CHROM == "X")] <- 23
data$CHROM <- as.numeric(data$CHROM)


data <- data[order(data$CHROM, data$POS),]
dim(data)

#summary(factor(data$FILTER))

## Get merged gene annotation
geneName <- ""
for(i in 1:nrow(data)){
  geneName[i] <- data[i, (min(which(!is.na(data[i, 5:(4 + NUMSAMPLES)])))+4)]
}

## Get existence matrix
varData <- data[,5:(4 + NUMSAMPLES)]

for(i in 1:nrow(varData)){
  varData[i,][which(!is.na(varData[i,]))] <- 1
  varData[i,][which(is.na(varData[i,]))] <- 0
}

colnames(varData) <- excel_sheets(FILENAME)
#View(varData)

varData <- data.frame(CHROM = data$CHROM,
                      POS = data$POS, 
                      varData)

#write.table(varData, "SupplementaryTableS1.VariantPatterns.txt", sep = "\t", quote=F, row.names = F)


## Drawing heatmap for each variant position with phenotype (Figure 5C)
visData <- varData[,3:13]
visData <- sapply(visData, as.numeric)

rownames(visData) <- data$ID



## Drawing Figure 5C
pheMetaData <- data.frame(fread("PheatmapMetadata.txt", sep = "\t"))
colnames(visData) <- pheMetaData$New           


annoPheatmap <- pheMetaData[,-c(1:2)]

rownames(annoPheatmap) <- pheMetaData$New
annoPheatmap$Subject <- factor(annoPheatmap$Subject)
annoPheatmap$Pathology <- factor(annoPheatmap$Pathology)
annoPheatmap$Malignancies <- factor(annoPheatmap$Malignancies)


colManual <- c("white", "royalblue1")

pheatmap(visData, clustering_distance_cols = "correlation", border_color = "black", clustering_method = "single",
         cluster_rows = F, show_rownames = F, col = colManual, annotation_col = annoPheatmap)

pheatmap(visData, clustering_distance_cols = "correlation", border_color = "black", clustering_method = "single",
         cluster_rows = F, show_rownames = F, col = colManual, annotation_col = annoPheatmap,
         filename = "test.pdf")

pheatmap(visData, clustering_distance_cols = "correlation", filename = "Figure5A.pdf")


## Draiwng Figure 5A
pca <- prcomp(t(visData))
pca

pcaData <- as.data.frame(pca$x[, 1:2]) # extract first two PCs
pcaData <- cbind(pcaData,
                 Subject = factor(pheMetaData$Subject),
                 Type = factor(pheMetaData$Type)) # add species to df

ggplot(pcaData) +
  aes(PC1, PC2, color = Subject, shape = Type) + # define plot area
  geom_point(size = 6) + # adding data points
  coord_fixed() +
  theme_classic2() +
  scale_color_d3()
  #theme(legend.position="none")

ggsave("Figure5A.pdf", width = 8, height = 4)


## Drawing Figure 5B

col4 <- colorRampPalette(c("grey", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
col <- c(rep("white", 48), col4(50))


M <- cor(visData)

colnames(M) <- paste0("S", 1:11)
rownames(M) <- paste0("S", 1:11)

pdf("Figure5B_upper.pdf", width = 6, height = 6)
corrplot(M, is.corr = FALSE, method = "ellipse",
         order = 'hclust', addrect = 4, col.lim = c(min(M), max(M)), tl.pos = 'd',
         col = col, tl.col = "black",
         type = 'upper')
dev.off()
#addCoef.col = 'black'

pdf("Figure5B_lower.pdf", width = 6, height = 6)
corrplot(M, is.corr = FALSE, method = 'number',
         order = 'hclust', addrect = 4, col.lim = c(min(M), max(M)), addCoef.col = 'black', tl.pos = 'd',
         col = col, tl.col = "black",
         type = 'lower') 
dev.off()


