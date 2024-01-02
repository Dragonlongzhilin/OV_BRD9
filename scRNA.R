#' @description: process the scRNA with CD45 isolate in OV-Mouse

library(Seurat)
library(patchwork)
library(clustree)
library(tidyverse)
library(plyr)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
set.seed(101)
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2)
set.resolutions <- seq(0.2, 1.2, by = 0.1)
setwd("/data/scRNA/Analysis")

FunctionScDblFinder.RNA <- function(Seurat.object, clusters = NULL, samples = NULL, nfeatures = 1000, dims = 20, dbr = NULL, n.core = 10, seed = 12345){
    set.seed(seed)
    require(scDblFinder)
    require(BiocParallel)
    sce <- as.SingleCellExperiment(Seurat.object)
    sce <- scDblFinder(sce, clusters = clusters, samples = samples, nfeatures = nfeatures, dims = dims, BPPARAM = MulticoreParam(n.core))
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.class, "scDblFinder.class")
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.score, "scDblFinder.score")
    return(Seurat.object)
}

#### ---- 1.load data ---- ####
samples <- c("KO_ovX3", "NT_ovX3")
scRNA.list <- sapply(samples, function(x){
  scRNA.data <- Read10X_h5(filename = paste0("/data/CellRangerResults/", x, "/outs/CellBender/cellbender_output_file_filtered_seurat.h5"))
  scRNA.data <- CreateSeuratObject(counts = scRNA.data, project = x, min.cells = 3, min.features = 200)
  scRNA.data[["mt_ratio_RNA"]] <- PercentageFeatureSet(scRNA.data, pattern = "^mt-")
  scRNA.data[["rp_ratio_RNA"]] <- PercentageFeatureSet(scRNA.data, pattern = "^Rpl|^Rps")
  scRNA.data <- subset(scRNA.data, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & mt_ratio_RNA < 10) # remove low-quality cells
  scRNA.data <- FunctionScDblFinder.RNA(Seurat.object = scRNA.data, nfeatures = 2000, dims = 30, n.core = 16)
  scRNA.data <- subset(scRNA.data, subset = scDblFinder.class == "singlet") # remove doublets
  return(scRNA.data)
})

scRNA <- merge(scRNA.list[[1]], y = scRNA.list[2:length(scRNA.list)],
               add.cell.ids = c("OV_KO", "OV_NT"),
               project = "scRNA")
orig.ident <- mapvalues(scRNA$orig.ident, 
                        from = samples, 
                        to = c("OV_KO", "OV_NT"))
scRNA <- AddMetaData(object = scRNA, metadata = orig.ident, col.name = "orig.ident")
scRNA.pro <- subset(scRNA, subset = nFeature_RNA < 8000 & nCount_RNA < 50000) # 9590 cells
scRNA.pro$orig.ident <- factor(scRNA.pro$orig.ident, levels = c("OV_NT", "OV_KO"))

#### remove the cells not exress Ptprc
Ptprc <- scRNA.pro@assays$RNA@counts["Ptprc",]
Ptprc.label <- rep("YES", length(Ptprc))
Ptprc.label[which(Ptprc == 0)] <- "NO"
scRNA.pro$Ptprc_exp <- Ptprc.label
scRNA.pro <- subset(scRNA.pro, subset = Ptprc_exp == "YES") # 9179 cells

color.sample <- c("#3C5488", "#E64B35")
names(color.sample) <- c("OV_NT", "OV_KO")

VlnPlot(
  object = scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"),
  ncol = 4,
  cols = color.sample,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# log10
VlnPlot(
  object = scRNA.pro,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 3,
  cols = color.sample,
  pt.size = 0,
  group.by = 'orig.ident'
) & xlab("") & yscale("log10", .format = TRUE) & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#### ---- All cells ---- ####
scRNA.pro <- NormalizeData(scRNA.pro)
scRNA.pro <- SCTransform(scRNA.pro, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE)
scRNA.pro <- RunPCA(scRNA.pro, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA.pro, ndims = 50)
scRNA.pro <- RunUMAP(scRNA.pro, reduction = "pca", dims = 1:30, verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>% FindClusters(resolution = set.resolutions, verbose = FALSE)

DimPlot(scRNA.pro, reduction = "umap", group.by = "orig.ident", cols = color.sample) + theme(legend.position = "top")
scRNA.pro$seurat_clusters <- scRNA.pro$SCT_snn_res.0.2
DimPlot(scRNA.pro, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, label.size = 5) + NoLegend()
DotPlot(scRNA.pro, features = features$Gene, assay = "SCT", group.by = "seurat_clusters", cols = c("#1e90ff", "#ff5a36"), dot.scale = 4) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 11), legend.text = element_text(size = 10), axis.text.y = element_text(size = 11))


#### differential analysis
scRNA.pro <- PrepSCTFindMarkers(scRNA.pro)
Idents(scRNA.pro) <- scRNA.pro$seurat_clusters
DEG.cluster <- FindAllMarkers(scRNA.pro, only.pos = TRUE,
                              group.by = "seurat_clusters",
                              assay = "SCT", min.pct = 0.25)
DEG.cluster.sig <- DEG.cluster[which(DEG.cluster$p_val_adj < 0.05 & DEG.cluster$avg_log2FC > 0.25),]
saveRDS(DEG.cluster.sig, file = "DEG.cluster.sig.rds")
write.table(DEG.cluster.sig, file = "DEG.cluster.sig.txt", quote = F, sep = "\t", row.names = F)

#### assign cell type
cluster.label <- as.character(scRNA.pro$seurat_clusters)
cluster.label <- gsub("^12$", "Mast cell", cluster.label)
cluster.label <- gsub("^13$", "Basophil", cluster.label)
cluster.label <- gsub("^11$", "pDC", cluster.label)
cluster.label <- gsub("^5$", "cDC", cluster.label)
cluster.label <- gsub("^10$", "B cell", cluster.label)
cluster.label <- gsub("^1$", "T cell", cluster.label)
cluster.label <- gsub("^4$", "NK cell", cluster.label)
cluster.label <- gsub("^3$", "Microglia", cluster.label)
cluster.label[which(cluster.label %in% c(0,2,6,8,7,9))] <- "Macrophage"
scRNA.pro <- AddMetaData(scRNA.pro, cluster.label, "cellType")
cellType.order <- c("Macrophage", "Microglia", "cDC", "pDC", "Mast cell", "Basophil", "T cell", "NK cell", "B cell")
scRNA.pro@meta.data$cellType <- factor(scRNA.pro@meta.data$cellType, levels = cellType.order)
saveRDS(scRNA.pro, file= "scRNA.rds")

color.cellType <- c("#E64B35", "#F39B7F", "#00A087", "#91D1C2", "#8491B4", "#7E6148", "#EE7B19", "#B8A00B", "#3C5488")
names(color.cellType) <- levels(scRNA.pro$cellType)
DimPlot(scRNA.pro, reduction = "umap", group.by = "orig.ident", cols = color.sample) + theme(legend.position = "top")
DimPlot(scRNA.pro, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, label.size = 5) + NoLegend()
DimPlot(scRNA.pro, reduction = "umap", cols = color.cellType, group.by = "cellType", label = T, repel = T, label.size = 5) + NoLegend()

cellRatio <- as.data.frame(table(scRNA.pro@meta.data[,c("orig.ident", "cellType")]))
cellNumber <- as.data.frame(table(scRNA.pro@meta.data[,c("orig.ident")]))
cellNumber <- rep(cellNumber$Freq, times = length(levels(scRNA.pro$cellType)))
cellRatio$cellNumber <- cellNumber
cellRatio$ratio <- round(cellRatio$Freq / cellRatio$cellNumber * 100, 2)
# KO/NT
OV.change <- cellRatio$ratio[seq(2, nrow(cellRatio), by = 2)] / cellRatio$ratio[seq(1, nrow(cellRatio), by = 2)]
OV.change <- data.frame(cellType = cellRatio$cellType[seq(2, nrow(cellRatio), by = 2)], change = OV.change)
OV.change$type <- "Increase"
OV.change$type[which(OV.change$change < 1)] <- "Decrease"
OV.change$type <- factor(OV.change$type, levels = c("Increase", "Decrease"))
OV.change$change <- log2(OV.change$change)
# remove cell type with cell number less than 50 cells
OV.change <- OV.change[which(!(OV.change$cellType %in% c("Mast cell", "Basophil"))),]
OV.change$cellType <- droplevels(OV.change$cellType)
ggbarplot(OV.change, x = "cellType", y = "change", fill = "type", color = "type", sort.val = "desc", palette = c("#E64B35", "#214DA9"), position = position_dodge(1), ylab = "Log2(fold change (KO vs NT))", xlab = "") + rotate_x_text(90)

#### ---- recluster T cell ---- ####
#### cluster (808 cells)
T.cellTypes <- levels(scRNA.pro$cellType)[grep("T cell", levels(scRNA.pro$cellType))]
scRNA.T <- subset(scRNA.pro, subset = cellType %in% T.cellTypes)
scRNA.T$cellType <- droplevels(scRNA.T$cellType)
DefaultAssay(scRNA.T) <- "RNA"
scRNA.T <- SCTransform(scRNA.T, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE)
scRNA.T <- RunPCA(scRNA.T, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA.T, ndims = 50)
scRNA.T <- RunUMAP(scRNA.T, reduction = "pca", dims = 1:25, verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:25, verbose = FALSE) %>% FindClusters(resolution = set.resolutions, verbose = FALSE)
DimPlot(scRNA.T, reduction = "umap", group.by = "orig.ident", cols = color.sample)
scRNA.T$seurat_clusters <- scRNA.T$SCT_snn_res.0.3
DimPlot(scRNA.T, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, label.size = 5) + NoLegend()

#### differential analysis
scRNA.T <- PrepSCTFindMarkers(scRNA.T)
Idents(scRNA.T) <- scRNA.T$seurat_clusters
DEG.Tcell <- FindAllMarkers(scRNA.T, only.pos = TRUE, group.by = "seurat_clusters", assay = "SCT", min.pct = 0.25)
DEG.Tcell.sig <- DEG.Tcell[which(DEG.Tcell$p_val_adj < 0.05 & DEG.Tcell$avg_log2FC > 0.25),]
saveRDS(DEG.Tcell.sig, file = "DEG.Tcell.sig.rds")
write.csv(DEG.Tcell.sig, file = "DEG.Tcell.sig.csv", quote = F, row.names = F)

#### cell type marker
T.markers <- c("Cd3d", "Cd3e", "Cd8a", "Cd4", "Foxp3", "Tcf7", "Lef1", "H2-Aa", "Cd86", "Il23r", "Il17a")
MySeuratWrappers::VlnPlot(scRNA.T, group.by = "seurat_clusters", features = T.markers, direction = "horizontal", stacked=T, pt.size=0, x.lab = "", y.lab = "")

#### assign the cell type
cluster.label <- as.character(scRNA.T$seurat_clusters)
cluster.label <- gsub("^0$", "CD8+ T-C1", cluster.label)
cluster.label <- gsub("^1$", "CD4+ T", cluster.label)
cluster.label <- gsub("^2$", "CD8+ T-C2", cluster.label)
cluster.label <- gsub("^3$", "Treg", cluster.label)
cluster.label <- gsub("^4$", "CD8+ T-C3", cluster.label)
cluster.label <- gsub("^5$", "Unknown", cluster.label)
cluster.label <- gsub("^6$", "Th17", cluster.label)
scRNA.T <- AddMetaData(scRNA.T, cluster.label, "cellType")
cellType.order <- c("CD8+ T-C1", "CD8+ T-C2", "CD8+ T-C3", "CD4+ T", "Treg", "Th17", "Unknown")
scRNA.T@meta.data$cellType <- factor(scRNA.T@meta.data$cellType, levels = cellType.order)
saveRDS(scRNA.T, "scRNA.T.rds")

color.T <- pal_npg("nrc")(length(unique(scRNA.T$cellType)))
names(color.T ) <- levels(scRNA.T$cellType)
DimPlot(scRNA.T, reduction = "umap", group.by = "orig.ident", cols = color.sample) + theme(legend.position = "top")
DimPlot(scRNA.T, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, label.size = 5) + NoLegend()
DimPlot(scRNA.T, reduction = "umap", group.by = "cellType", cols = color.T, label = T, repel = T, label.size = 5) + NoLegend()

cellRatio <- as.data.frame(table(scRNA.T@meta.data[,c("orig.ident", "cellType")]))
cellNumber <- as.data.frame(table(scRNA.T@meta.data[,c("orig.ident")]))
cellNumber <- rep(cellNumber$Freq, times = length(levels(scRNA.T$cellType)))
cellRatio$cellNumber <- cellNumber
cellRatio$ratio <- round(cellRatio$Freq / cellRatio$cellNumber * 100, 2)
OV.change <- cellRatio$ratio[seq(2, nrow(cellRatio), by = 2)] / cellRatio$ratio[seq(1, nrow(cellRatio), by = 2)]
OV.change <- data.frame(cellType = cellRatio$cellType[seq(2, nrow(cellRatio), by = 2)], change = OV.change)
OV.change$type <- "Increase"
OV.change$type[which(OV.change$change < 1)] <- "Decrease"
OV.change$type <- factor(OV.change$type, levels = c("Increase", "Decrease"))
OV.change$change <- log2(OV.change$change)
ggbarplot(OV.change, x = "cellType", y = "change", fill = "type", color = "type", sort.val = "desc", palette = c("#E64B35", "#214DA9"), position = position_dodge(1), ylab = "Log2(fold change (KO vs NT))", xlab = "") + rotate_x_text(90)

#### ---- Macrophage cell ---- ####
#### cluster (6897 cells)
Macrophage.cellTypes <- levels(scRNA.pro$cellType)[grep("Macrophage", levels(scRNA.pro$cellType))]
scRNA.Macrophage <- subset(scRNA.pro, subset = cellType %in% Macrophage.cellTypes)
scRNA.Macrophage@meta.data <- scRNA.Macrophage@meta.data[,-grep("SCT_snn_res", colnames(scRNA.Macrophage@meta.data))]
scRNA.Macrophage$cellType <- droplevels(scRNA.Macrophage$cellType)
DefaultAssay(scRNA.Macrophage) <- "RNA"
scRNA.Macrophage <- SCTransform(scRNA.Macrophage, vst.flavor = "v2", vars.to.regress = "nCount_RNA", verbose = FALSE)
scRNA.Macrophage <- RunPCA(scRNA.Macrophage, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA.Macrophage, ndims = 50)
scRNA.Macrophage <- RunUMAP(scRNA.Macrophage, reduction = "pca", dims = 1:30, verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>% FindClusters(resolution = seq(0.1, 0.6, by = 0.1), verbose = FALSE)
DimPlot(scRNA.Macrophage, reduction = "umap", group.by = "orig.ident", cols = color.sample)
DimPlot(scRNA.Macrophage, reduction = "umap", group.by = "SCT_snn_res.0.4", label = T, repel = T, label.size = 5) + NoLegend()
scRNA.Macrophage$seurat_clusters <- scRNA.Macrophage$SCT_snn_res.0.4
scRNA.Macrophage$cellType <- factor(paste0("C", as.numeric(scRNA.Macrophage$seurat_clusters)), levels = paste0("C", 1:max(as.numeric(scRNA.Macrophage$seurat_clusters))))
saveRDS(scRNA.Macrophage, "scRNA.Macrophage.rds")

#### differential analysis
scRNA.Macrophage <- PrepSCTFindMarkers(scRNA.Macrophage)
Idents(scRNA.Macrophage) <- scRNA.Macrophage$cellType
DEG.Macrophage <- FindAllMarkers(scRNA.Macrophage, only.pos = TRUE, group.by = "cellType", assay = "SCT", min.pct = 0.25)
DEG.Macrophage.sig <- DEG.Macrophage[which(DEG.Macrophage$p_val_adj < 0.05 & DEG.Macrophage$avg_log2FC > 0.25),]
saveRDS(DEG.Macrophage.sig, file = "DEG.Macrophage.sig.rds")
write.csv(DEG.Macrophage.sig, file = "DEG.Macrophage.sig.csv", quote = F, row.names = F)

color.M <- pal_npg("nrc")(length(unique(scRNA.Macrophage$cellType)))
color.M[9] <- "#EE7B19"
names(color.M) <- levels(scRNA.Macrophage$cellType)
DimPlot(scRNA.Macrophage, reduction = "umap", group.by = "orig.ident", cols = color.sample) + theme(legend.position = "top")
DimPlot(scRNA.Macrophage, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, label.size = 5) + NoLegend()
DimPlot(scRNA.Macrophage, reduction = "umap", group.by = "cellType", cols = color.M, label = T, repel = T) + NoLegend()

#### cell ratio plot
cellRatio <- as.data.frame(table(scRNA.Macrophage@meta.data[,c("orig.ident", "cellType")]))
cellNumber <- as.data.frame(table(scRNA.Macrophage@meta.data[,c("orig.ident")]))
cellNumber <- rep(cellNumber$Freq, times = length(levels(scRNA.Macrophage$cellType)))
cellRatio$cellNumber <- cellNumber
cellRatio$ratio <- round(cellRatio$Freq / cellRatio$cellNumber * 100, 2)
OV.change <- cellRatio$ratio[seq(2, nrow(cellRatio), by = 2)] / cellRatio$ratio[seq(1, nrow(cellRatio), by = 2)]
OV.change <- data.frame(cellType = cellRatio$cellType[seq(2, nrow(cellRatio), by = 2)], change = OV.change)
OV.change$type <- "Increase"
OV.change$type[which(OV.change$change < 1)] <- "Decrease"
OV.change$type <- factor(OV.change$type, levels = c("Increase", "Decrease"))
OV.change$change <- log2(OV.change$change)
ggbarplot(OV.change, x = "cellType", y = "change", fill = "type", color = "type", sort.val = "desc", palette = c("#E64B35", "#214DA9"), position = position_dodge(1), ylab = "Log2(fold change (KO vs NT))", xlab = "") + rotate_x_text(90)
