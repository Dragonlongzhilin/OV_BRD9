#' @description: analysis RNA-seq

setwd("/data/RNA-seq/Analysis")
#### ---- 1.load data ---- ####
countMatrix <- readRDS("countMatrix.rds")

#### ---- differential analysis ---- ####
sampleNames <- colnames(countMatrix)
phenotype.info <- data.frame(SampleName = c("NTm-1", "NTm-2", "B29m-1", "B29m-2", "B30m-1", "B30m-2", "NT3-1", "NT3-2", "B293-1", "B293-2", "B303-1", "B303-2"),
							 condition = c("Mock-NT", "Mock-NT", "Mock-KO", "Mock-KO", "Mock-KO", "Mock-KO", "OV-NT", "OV-NT", "OV-KO", "OV-KO", "OV-KO", "OV-KO"))
condition <- plyr::mapvalues(sampleNames, phenotype.info$SampleName, phenotype.info$condition)
coldata <- data.frame(condition = condition)
rownames(coldata) <- sampleNames

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=coldata, design=~ condition)
dds <- DESeq(dds)
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("condition"))

Mock_NT.vs.KO <- results(dds, contrast=c("condition", "Mock-NT", "Mock-KO"))
Mock_NT.vs.KO <- as.data.frame(Mock_NT.vs.KO)
Mock_NT.vs.KO <- na.omit(Mock_NT.vs.KO)
Mock_NT.vs.KO$gene <- rownames(Mock_NT.vs.KO)

OV_NT.vs.KO <- results(dds, contrast=c("condition", "OV-NT", "OV-KO"))
OV_NT.vs.KO <- as.data.frame(OV_NT.vs.KO)
OV_NT.vs.KO <- na.omit(OV_NT.vs.KO)
OV_NT.vs.KO$gene <- rownames(OV_NT.vs.KO)

KO_Mock.vs.OV <- results(dds, contrast=c("condition", "Mock-KO", "OV-KO"))
KO_Mock.vs.OV <- as.data.frame(KO_Mock.vs.OV)
KO_Mock.vs.OV <- na.omit(KO_Mock.vs.OV)
KO_Mock.vs.OV$gene <- rownames(KO_Mock.vs.OV)

NT_Mock.vs.OV <- results(dds, contrast=c("condition", "Mock-NT", "OV-NT"))
NT_Mock.vs.OV <- as.data.frame(NT_Mock.vs.OV)
NT_Mock.vs.OV <- na.omit(NT_Mock.vs.OV)
NT_Mock.vs.OV$gene <- rownames(NT_Mock.vs.OV)

diff.merge.result <- list(Mock_NT.vs.KO = Mock_NT.vs.KO, 
						  OV_NT.vs.KO = OV_NT.vs.KO,
						  KO_Mock.vs.OV = KO_Mock.vs.OV,
						  NT_Mock.vs.OV = NT_Mock.vs.OV)
saveRDS(diff.merge.result, file = "diff.merge.result.rds")
writexl::write_xlsx(diff.merge.result, path = "diff.merge.result.xlsx")

# significant change: abs(log2(FC)) > 1 & FDR < 0.05
diff.merge.result.sig <- lapply(diff.merge.result, function(x){
	x.up <- x[which(x$log2FoldChange > 1 & x$padj < 0.05),]
	x.down <- x[which(x$log2FoldChange < -1 & x$padj < 0.05),]
	return(list(Up = data.frame(gene = rownames(x.up)), Down = data.frame(gene = rownames(x.down))))
})
diff.merge.result.sig <- unlist(diff.merge.result.sig, recursive = F)
writexl::write_xlsx(diff.merge.result.sig, path = "diff.merge.result.sig.xlsx")

#### Figure heatmap
library(ComplexHeatmap)
library(circlize)
log2TPM <- readRDS("log2TPM.rds")
heatmapgenes <- read.table("/data/heatmapGene.txt", header = T, stringsAsFactors = F, sep = "\t")

# order columns
log2TPM <- log2TPM[,c(11,12,3:4,7:8,9,10,1:2,5:6)]
coldata <- coldata[c(11,12,3:4,7:8,9,10,1:2,5:6),,drop = F]

col_split <- gsub("-.*", "", coldata$condition)
column_ha <- HeatmapAnnotation(condition = coldata$condition, annotation_legend_param = list(legend_direction = "horizontal"), col = list(condition = c("OV-KO" = "#E64B35", "OV-NT" = "#00A087", "Mock-KO" = "#E18727", "Mock-NT" = "#3C5488")))

idx <- match(heatmapgenes$Gene, rownames(log2TPM))
row_split <- heatmapgenes$Type
plot.data.TPM <- t(scale(t(log2TPM[idx,])))
colnames(plot.data.TPM) <- c("NT-1", "NT-2", "KO-1-1", "KO-1-2", "KO-2-1", "KO-2-2", "NT-1", "NT-2", "KO-1-1", "KO-1-2", "KO-2-1", "KO-2-2")
Heatmap(plot.data.TPM, name = "Scaled expression", cluster_columns = F, column_split = col_split, row_split = row_split, top_annotation = column_ha, width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
