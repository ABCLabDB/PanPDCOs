#################################
## scRNA-seq for PDAC organoid ##
## Dimensional Reudction       ##
## Minseok Seo                 ##
#################################

## Set Repositories
setRepositories(ind = 1:7)

## Install packages
#BiocManager::install("harmony", version = "3.14")
#devtools::install_github("thomasp85/patchwork")
#install.packages("scCustomize")

## Load library
library(Seurat)
library(patchwork)
library(harmony)
library(cowplot)
library(data.table)
library(cluster)
library(factoextra)
library(future)
library(ggpubr)
library(dplyr)
library(writexl)
library(preprocessCore)
library(scCustomize)
library(wesanderson)

## Set Work & Data dir.
WORK_DIR <- ""
DATA_DIR <- "./0.Data/1.scRNASeq"

## Load raw count data
setwd(DATA_DIR)

firstcount <- t(fread("PO1_SampleTag01_hs_128B_RSEC_MolsPerCell.csv"))
colnames(firstcount) <- cellID <- firstcount[1,]
firstcount <- firstcount[-1,]

firstObject <- CreateSeuratObject(counts = firstcount , project = "Tag1")

fileList <- list.files(path = "./", pattern = "MolsPerCell.csv")
projectID <- c("Tag2", "Tag3","Tag4", "Tag5","Tag7","Tag8","Tag9")

object.list <- vector(mode = "list", length = length(fileList)-1)

for (i in 2:length(fileList)) {
  
  temp <- t(fread(fileList[i]))
  colnames(temp) <- cellID <- temp[1,]
  temp <- temp[-1,]
  
  object.list[[(i-1)]] <- CreateSeuratObject(counts = temp, project =  projectID[(i-1)])
  print(i)
}

object.combined <- merge(firstObject, y=object.list, add.cell.ids = c("Tag1", projectID), project="AllObject")


## Check loaded data
object.combined
head(colnames(object.combined))
tail(colnames(object.combined))
table(object.combined$orig.ident)
unique(sapply(X = strsplit(colnames(object.combined), split = "_"), FUN = "[", 1))


## Loaded meta information
meta_data <- data.frame(fread("Metadata.txt"))
object.combined$library <- object.combined$orig.ident
object.combined$library <- plyr::mapvalues(object.combined$library, meta_data$Sample.Tag, c("PO1","PO1","PO1","PO4","PO4","PO7","PO7","PO7"))
object.combined$library <- as.factor(object.combined$library)
levels(object.combined$library)


## Add Sampling methods for meta-data
object.combined$sampleMethods <- object.combined$orig.ident
object.combined$sampleMethods <- plyr::mapvalues(object.combined$sampleMethods, meta_data$Sample.Tag, meta_data$Sample.Methods)
object.combined$sampleMethods <- as.factor(object.combined$sampleMethods)

## Add Diagnosis for meta-data
object.combined$diagnosis <- object.combined$orig.ident
object.combined$diagnosis <- plyr::mapvalues(object.combined$diagnosis, meta_data$Sample.Tag, meta_data$Etc)
object.combined$diagnosis <- as.factor(object.combined$diagnosis)

## Dx for meta-data
object.combined$DX <- object.combined$orig.ident
object.combined$DX <- plyr::mapvalues(object.combined$DX, meta_data$Sample.Tag, meta_data$DX)
object.combined$DX <- as.factor(object.combined$DX)

## Interaction
object.combined$group <- object.combined$orig.ident
object.combined$group <- plyr::mapvalues(object.combined$group, meta_data$Sample.Tag, as.character(factor(meta_data$DX):factor(meta_data$Etc))
)
object.combined$group <- as.factor(object.combined$group)


##Pre-processing (before batch correction)
object.combined[["percent.mt"]] <- PercentageFeatureSet(object = object.combined, pattern = "^MT-")
object.combined <- subset(object.combined, subset = nFeature_RNA > 400 & percent.mt < 40)
object.combined <- NormalizeData(object.combined) %>% 
  FindVariableFeatures(selection.method='vst', nfeatures=2000) %>% 
  ScaleData() %>% 
  RunPCA()


## HARMONY correction
object.combined <- object.combined %>% 
  RunHarmony("library", plot_convergence=TRUE)

harmony_embeddings <- Embeddings(object.combined, 'harmony')

## For DOWNSTREAM ANALYSIS
object.Umap <- object.combined
object.tSNE <- object.combined


## UMAP clustering
object.Umap <- object.Umap %>% 
  RunUMAP(reduction='harmony', dims=1:10) %>% 
  FindNeighbors(reduction='harmony', dims=1:10) %>% 
  FindClusters(resolution=0.1) %>% 
  identity()

distance_matrix <- dist(Embeddings(object.Umap[['umap']])[, 1:2])
clusters <- object.Umap@active.ident
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
object.Umap@meta.data$silhouette_score <- silhouette[,3]
fviz_silhouette(silhouette)

## Drawing Figures
setwd("/disk2/bisms/2.scRNASeqPDAC_sms/3.Figures")

pdf("Figure4A.UMAPCluster_6group_noLgend.pdf", height = 4, width = 4)
my_cols = DiscretePalette_scCustomize(num_colors = 6, palette = "varibow")
DimPlot(object.Umap, reduction = 'umap', label = TRUE, cols=alpha(my_cols,0.3), pt.size = .8) + NoLegend()
dev.off()


pdf("Figure4D.legend.pdf", height = 4, width = 4)
my_cols = JCO_Four()
DimPlot(object.Umap, reduction = 'umap', group.by = 'diagnosis', cols=alpha(my_cols,0.3), pt.size = .8)
dev.off()


pdf("Figure4D.DiagnosisUMAP_noLegend.pdf", height = 4, width = 4)
my_cols = JCO_Four()
DimPlot(object.Umap, reduction = 'umap', group.by = 'diagnosis', cols=alpha(my_cols,0.3), pt.size = .8) + NoLegend()
dev.off()


### Cluster biomarkers (finding differentially expressed features in each cluster)
object.markers <- FindAllMarkers(object.Umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
toppercluster <- object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

setwd(WORK_DIR)
#write_xlsx(toppercluster, path = "./Harmony_markers.xlsx")


## Figure 4B. Drawing Heatmap
object.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#DoHeatmap(object.Umap, features = top5$gene) + scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))
DoHeatmap(object.Umap, features = top5$gene, lines.width = 2) + scale_fill_gradientn(colours = myPalette)
ggplot2::ggsave(filename = "Figure4B.Legend.pdf", height = 6, width = 4)

DoHeatmap(object.Umap, features = top5$gene, lines.width = 2) + scale_fill_gradientn(colours = myPalette) + NoLegend()
ggplot2::ggsave(filename = "Figure4B.NoLegend.pdf", height = 6, width = 4)


## Figure 4C. Drawing gene plots
setwd("/disk2/bisms/2.scRNASeqPDAC_sms/3.Figures")

my_cols = DiscretePalette_scCustomize(num_colors = 6, palette = "varibow")

VlnPlot(object.Umap, features = c("TFF1","KRT19","S100A10",
                                  "LGALS4","PERP","CLDN4"), ncol = 2,
        cols=alpha(my_cols,0.3), pt.size = .8) + NoLegend()

#ggplot2::ggsave(filename = "Figure4C.DuctaclCellMarkers.pdf")


markerName <- "TFF1"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)



markerName <- "KRT19"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)


markerName <- "S100A10"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)



markerName <- "LGALS4"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)




markerName <- "PERP"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)





markerName <- "CLDN4"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Cluster", y="Expression", palette = "d3") %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.05, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Cluster", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") +
  NoLegend()

ggsave(paste0("Figure4C.", markerName, ".pdf"), width = 5, height = 3)



## Get Avg Expression Matrix from whole cells
avgExpMatrix <- AverageExpression(object.Umap)
avgExpMatrix <- avgExpMatrix$RNA
avgExpMatrix <- data.frame(WholeCellAvgExp = rowSums(avgExpMatrix))
avgExpMatrix <- data.frame(geneSymbol = rownames(avgExpMatrix),
                           avgExpMatrix)


## For subset of cells from UMAP outcome
object.Umap.PathNeg <- subset(object.Umap, subset = diagnosis == "Pathology (-)")
object.Umap.PathPos <- subset(object.Umap, subset = diagnosis == "Pathology (+)")
object.Umap.Suspicious <- subset(object.Umap, subset = diagnosis == "Suspicious")
object.Umap.GBC <- subset(object.Umap, subset = DX == "GBC")
object.Umap.PDAC <- subset(object.Umap, subset = DX == "PDAC")

object.Umap.GBC <- subset(object.Umap, subset = DX == "GBC")
object.Umap.PDAC <- subset(object.Umap, subset = DX == "PDAC")

object.Umap.C0 <- subset(object.Umap, idents = '0')
object.Umap.C1 <- subset(object.Umap, idents = '1')
object.Umap.C2 <- subset(object.Umap, idents = '2')
object.Umap.C3 <- subset(object.Umap, idents = '3')
object.Umap.C4 <- subset(object.Umap, idents = '4')
object.Umap.C5 <- subset(object.Umap, idents = '5')


## Drawing Figure 5

fig5List <- c("LYZ",
              "TFF1",
              "TFF2",
              "EEF1A1",
              "S100A6",
              "TPT1",
              "PTMA",
              "LCN2",
              "ANXA2",
              "AGR2",
              "ACTG1",
              "CFL1",
              "B2M",
              "SERPINA1",
              "PPIA",
              "CLDN2",
              "TMSB4X",
              "KRT8",
              "MALAT1",
              "FTH1")

## Drawing Figure 5A (Drawing malignant ductal cells from whole cells)
FeaturePlot(object.Umap, features = fig5List, ncol = 5, keep.scale = "all") & NoAxes() & NoLegend() +
  theme( plot.title = element_text(face = "italic", size = 10))
ggsave("Figure5A.pdf", width=8, height=4)

FeaturePlot(object.Umap, features = fig5List, ncol = 5, keep.scale = "all") & NoAxes()
ggsave("Figure5A_legend.jpg", width=8, height=8)


## Drawing Figure 5B
FeaturePlot(object.Umap.PathNeg, features = fig5List, ncol = 5, keep.scale = "all", col = rev(viridis_plasma_dark_high)) & NoLegend() & NoAxes() +
  theme( plot.title = element_text(face = "italic", size = 10))

ggsave("Figure5B.pdf", width=8, height=4)


FeaturePlot(object.Umap.PathNeg, features = fig5List, ncol = 5, keep.scale = "all", col = rev(viridis_plasma_dark_high)) 
ggsave("Figure5B_legend.jpg", width=8, height=8)



## Figure 5C
markerName <- "LYZ"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Diagnosis = object.Umap@meta.data$diagnosis)

temp$Diagnosis <- factor(temp$Diagnosis, level = c("Pathology (-)", "Suspicious", "Pathology (+)"))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Diagnosis", y="Expression", palette = "aaas") %>% 
  ggadd(c("jitter"), color = "Diagnosis", alpha = 0.03, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Diagnosis", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5C.", markerName, ".pdf"), width = 6, height = 3)



## Figure 5D
markerName <- "TFF1"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Diagnosis = object.Umap@meta.data$diagnosis)
temp <- temp[-which(temp$Expression == 0),]


temp$Diagnosis <- factor(temp$Diagnosis, level = c("Pathology (-)", "Suspicious", "Pathology (+)"))

ggviolin(temp, x="Diagnosis", y="Expression", palette = "aaas") %>% 
  ggadd(c("jitter"), color = "Diagnosis", alpha = 0.03, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Diagnosis", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5D.", markerName, ".pdf"), width = 6, height = 3)




## Figure 5E
markerName <- "EEF1A1"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Diagnosis = object.Umap@meta.data$diagnosis)

temp$Diagnosis <- factor(temp$Diagnosis, level = c("Pathology (-)", "Suspicious", "Pathology (+)"))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Diagnosis", y="Expression", palette = "aaas") %>% 
  ggadd(c("jitter"), color = "Diagnosis", alpha = 0.03, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Diagnosis", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5E.", markerName, ".pdf"), width = 6, height = 3)


## Figure 5F
markerName <- "TPT1"
temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Diagnosis = object.Umap@meta.data$diagnosis)

temp$Diagnosis <- factor(temp$Diagnosis, level = c("Pathology (-)", "Suspicious", "Pathology (+)"))
temp <- temp[-which(temp$Expression == 0),]

ggviolin(temp, x="Diagnosis", y="Expression", palette = "aaas") %>% 
  ggadd(c("jitter"), color = "Diagnosis", alpha = 0.03, size = 3) %>% 
  ggadd(c("boxplot"), fill = "Diagnosis", alpha = 0.1) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5F.", markerName, ".pdf"), width = 6, height = 3)



##  Figure 5G

pdf("Figure5G.UMAPCluster.pdf", height = 3, width = 6)
my_cols = NavyAndOrange()
DimPlot(object.Umap, reduction='umap', group.by = 'DX', cols=alpha(my_cols, 0.3), pt.size = .8) + NoLegend() +  ggtitle(NULL)
dev.off()

pdf("Figure5G.UMAPCluster_legend.pdf", height = 3, width = 6)
my_cols = NavyAndOrange()
DimPlot(object.Umap, reduction='umap', group.by = 'DX', cols=alpha(my_cols, 0.3), pt.size = .8)
dev.off()


## Figure 5H
markerName <- "VCAN"
my_cols = DiscretePalette_scCustomize(num_colors = 6, palette = "varibow")

temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))

ggboxplot(temp, x="Cluster", y="Expression", palette = my_cols, outlier.shape = NA) %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.03, size = 3) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5H.", markerName, ".pdf"), width = 6, height = 3)


## Figure 5I
markerName <- "AQP3"

temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)

temp$Cluster <- factor(paste0("c", temp$Cluster))

ggboxplot(temp, x="Cluster", y="Expression", palette = my_cols, outlier.shape = NA) %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.03, size = 3) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5I.", markerName, ".pdf"), width = 6, height = 3)


## Figure 5J
markerName <- "FGF19"

temp <- data.frame(Expression = FetchData(object.Umap, vars = c(markerName))[,1],
                   Cluster = object.Umap@active.ident)
temp$Cluster <- factor(paste0("c", temp$Cluster))

ggboxplot(temp, x="Cluster", y="Expression", palette = my_cols, outlier.shape = NA) %>% 
  ggadd(c("jitter"), color = "Cluster", alpha = 0.03, size = 3) +
  xlab("") +
  ylab("Expression level") & NoLegend()

ggsave(paste0("Figure5J.", markerName, ".pdf"), width = 6, height = 3)

#



## FIN


