install.packages("SingleR")
install.packages("SingleCellExperiment", dependencies = TRUE)
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(hdf5r)
library(openxlsx)
library(readxl)
library(rio)
library(data.table)
library(magrittr)
library(utils)
library(SingleR)
library("SingleCellExperiment")

rewildingproject <- readRDS(file = "rewildingprojectobject.rds")

#Integrating Block 1 and Block 2

rewilding <- readRDS(file = "/Users/howardng/Documents/rewilding2.rds")
rewilding.combined <- readRDS(file = "/Users/howardng/Documents/CellTypexEnvironIdents_RewildingCombined.rds")

#Try combining the datasets now. Only make one UMAP and then split by treatment for analysis.
rewildingblocks.anchors <- FindIntegrationAnchors(object.list = list(rewilding, rewilding.combined), dims = 1:30)
rewildingproject <- IntegrateData(anchorset = rewildingblocks.anchors, dims = 1:30)
DefaultAssay(rewildingproject) <- "integrated"
rewildingproject <- NormalizeData(rewildingproject)
rewildingproject <- ScaleData(rewildingproject)
rewildingproject <- NormalizeData(rewildingproject, assay = "RNA")
rewildingproject <- ScaleData(rewildingproject, assay = "RNA")
rewildingproject <- NormalizeData(rewildingproject, assay = "ADT", normalization.method = "CLR")
rewildingproject <- ScaleData(rewildingproject, assay = "ADT")
rewildingproject <- NormalizeData(rewildingproject, assay = "HTO", normalization.method = "CLR")
rewildingproject <- ScaleData(rewildingproject, assay = "HTO")

saveRDS(rewildingproject, file = "/Users/howardng/documents/rewildingproject.rds")

table(rewildingproject$HTO_classification.global)
table(rewildingproject$HTO_maxID)
table(rewildingproject$stim)
View(rewildingproject@meta.data)

rewildingproject <- FindVariableFeatures(rewildingproject, selection.method = "vst", nfeatures = 5000)
top10RNA <- head(VariableFeatures(rewildingproject), 10)
View(top10RNA)
plot1 <- VariableFeaturePlot(rewildingproject)
plot2 <-LabelPoints(plot = plot1, points = top10RNA, repel = TRUE)
plot2

rewildingproject <- FindVariableFeatures(rewildingproject, assay = "ADT")
top10ADT <- head(VariableFeatures(rewildingproject$ADT), 10)
View(top10ADT)
plot1 <- VariableFeaturePlot(rewildingproject, assay = "ADT")
plot2 <-LabelPoints(plot = plot1, points = top10ADT, ynudge = 0,repel = TRUE)
plot2


####Clustering

#PCA
rewildingproject <- RunPCA(rewildingproject)
print(rewildingproject[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(rewildingproject, dims = 1:2, reduction = "pca")
DimPlot(rewildingproject, reduction = "pca")
DimHeatmap(rewildingproject, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(rewildingproject, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(rewildingproject, ndims = 40, reduction = "pca")


#cluster the cells
rewildingproject <- FindNeighbors(rewildingproject, dims = 1:20)
rewildingproject <- FindClusters(rewildingproject, resolution = 0.2)

#try .2 - .8 cluster resolutions to see what makes a better umap

#run non-linear dimensional reduction (UMAP/tSNE)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.3)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.4)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.5)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.6)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.7)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

rewildingproject <- FindClusters(rewildingproject, resolution = 0.8)
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")

#UMAPs
Idents(rewildingproject) <- "environmentgroup"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "infectionstatus"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "strain"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "Block"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "orig.ident"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")

saveRDS(rewildingproject, file = "rewildingprojectobject.rds")

rewildingproject <- readRDS(file = "/Users/howardng/documents/Rewilding Project/Block 1 + Block 2/objects/RPclustered.rds")

View(rewildingproject@meta.data)

DefaultAssay(rewildingproject) <- "RNA"

##Find Markers
Idents(rewildingproject) <- "seurat_clusters"
#RNA
RP.rna.markers <- FindAllMarkers(rewildingproject, only.pos = T)
RP.rna.markers <- RP.rna.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
DoHeatmap(rewildingproject, RP.rna.markers$gene, size = 3)+
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#6300ae','#F8766D','#FFCC99')), midpoint = 0, guide = "colourbar", aesthetics = "fill") 
RP.rna.markers <- RP.rna.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
export(RP.rna.markers, file = "/Users/howardng/documents/Rewilding Project/Block 1 + Block 2/RNAassay_markers_top25.xlsx")

#ADT
RP.adt.markers <- FindAllMarkers(rewildingproject, assay ="ADT", only.pos = T)
RP.adt.markers <- RP.adt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(rewildingproject, RP.adt.markers$gene, size = 3, assay = "ADT")+
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#6300ae','#F8766D','#FFCC99')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
export(RP.adt.markers, file = "/Users/howardng/documents/Rewilding Project/Block 1 + Block 2/Block1and2_ADTmarkers.xlsx")


Idents(rewildingproject) <- "orig.ident"
rewildingproject@meta.data$batch <- ifelse(rewildingproject@meta.data$orig.ident %in% c("rewild7.hashtag", "rewild8.hashtag", "rewild9.hashtag", "rewild10.hashtag",
                                                                                        "rewild11.hashtag", "rewild12.hashtag", "rewild13.hashtag", "rewild14.hashtag"), "Batch2", "Batch1")
options(ggrepel.max.overlaps = Inf)

View(rewildingproject@meta.data)

Idents(rewildingproject) <- "CellType"
rewildingprojectPCA <- sc_pca(rewildingproject, "CellType", "hash.ID", 
                              "batch")
ggsave("scPCA_batch.pdf", rewildingprojectPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/")

RPenvironmentPCA <- sc_pca(rewildingproject, "seurat_clusters", "hash.ID", "environmentgroup")
ggsave("RP_scPCA_environment.pdf", RPenvironmentPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/")

RPinfectionPCA <- sc_pca(rewildingproject, "seurat_clusters", "hash.ID", "infectionstatus")
ggsave("RP_scPCA_infection3.pdf", RPinfectionPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/")

RPstrainPCA <- sc_pca(rewildingproject, "seurat_clusters", "hash.ID", "strain")
ggsave("RP_scPCA_strain.pdf", RPstrainPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/")

Idents(rewildingproject)

Idents(rewildingproject) <- "infectionstatus"
rewildingproject <- RenameIdents(rewildingproject, "Not Infected" = "Uninfected")

RPinfectionstatusPCA <- sc_pca(rewildingproject, "seurat_clusters", "hash.ID", "infectionstatus")
ggsave("RP_scPCA_infectionstatus.pdf", RPinfectionstatusPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/")

#Separating Experimental Groups by Environment
rewilding@meta.data$environmentgroup <- ifelse(rewilding@meta.data$HTO_maxID %in% c("296", "293", "281", "290", "271", "286",
                                                                                    "299", "272", "285", "279", "306", "303", "313.323A",
                                                                                    "325", "321", "301", "312", "313.323B", "308", "319",
                                                                                    "102", "113", "105", "116", "119", "112", "121", "115",
                                                                                    "101", "107"), "Lab", "Rewilded")


Idents(rewilding)<- "environmentgroup"
View(rewilding@meta.data)



rewildingproject@meta.data$infectionstatus <- ifelse(rewildingproject@meta.data$hash.ID %in% c("434", "428", "445", "435", "478",
                                                                                                   "446", "427", "7010", "415", "402",
                                                                                                   "408", "412", "459", "450", "465",
                                                                                                   "403", "405", "410", "417", "422",
                                                                                                   "454", "456", "461", "468", "114", "104", "111", "123", "33",
                                                                                                   "320", "307", "304", "296", "293", "281",
                                                                                                   "290", "271", "306", "303", "313.323A", 
                                                                                                   "325", "321", "102", "113", "105", "116",
                                                                                                   "119", "275", "298", "291", "294", "287", "283",
                                                                                                   "280", "274", "311"), "Uninfected", "Infected")

rewildingproject@meta.data$infectionstatus <- factor(rewildingproject@meta.data$infectionstatus, levels=c("Uninfected", "Infected"))
is.factor(rewildingproject@meta.data$infectionstatus)
table(rewildingproject$infectionstatus)

View(rewildingproject@meta.data)

saveRDS(rewildingproject, file = "RewildingProject.rds")

#Using SingleR to identify clusters
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(SingleR)

head(rewildingproject)
RPsubset <- head(rewildingproject)
pred.RPsubset <- SingleR(test = RPsubset, ref = mpa.se, assay.type.test=1,
                     labels = mpa.se$label.main)
table(mpa.se$label.main)

imgen.se <- ImmGenData()
imgen.se
table(imgen.se$label.main)

write.table(as.matrix(GetAssay(object = rewildingproject, slot = "counts")), 
            '/Users/howardng/Documents/Rewilding Project/counts2.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

GetAssay(object = rewildingproject, slot = "counts")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
library("SingleCellExperiment")


RPmatrix <- write.table(rewildingproject@assays[["RNA"]]@counts, file='Counts.tsv', quote=FALSE, sep='\t', col.names = TRUE)
RPmatrix <- read.table(file = "Counts.tsv")
head(RPmatrix)
RPmatrixsub <- RPmatrix[1:10,1:10]
pred.RPmatrixsub <- SingleR(test = RPmatrixsub, ref = imgen.se, assay.type.test=1,
                         labels = imgen.se$label.main)
pred.RPmatrixsub
table(pred.RPmatrixsub$labels)


pred.RPmatrix <- SingleR(test = RPmatrix, ref = imgen.se, assay.type.test=1,
                         labels = imgen.se$label.main)
table(pred.RPmatrix$labels)

rewildingproject[["singleR main label"]] <- Idents(object = pred.RPmatrix)

Idents(rewilding.combined)<- "hash.ID"
rewilding.combined[["singleR main label"]] <- Idents(object = rewilding.combined)
Idents(rewilding.combined)<- 'sampleID'

predfine.RPmatrix <- SingleR(test = RPmatrix, ref = imgen.se, assay.type.test=1,
                         labels = imgen.se$label.fine)
table(predfine.RPmatrix$labels)
export(table(predfine.RPmatrix$labels), file = "SingleRfine.xlsx")
export(table(pred.RPmatrix$labels), file = "SingleRmain.xlsx")
export(pred.RPmatrix, file = "SingleRmainMatrix.xlsx")
export(predfine.RPmatrix, file = "SingleRfineMatrix.xlsx")
rownames(pred.RPmatrix)
export(rownames(pred.RPmatrix), file = "Rows_mainmatrix.xlsx")
export(rownames(predfine.RPmatrix), file = "Rows_finematrix.xlsx")



rewildingproject <- AddMetaData(object = rewildingproject, metadata= pred.RPmatrix$labels, col.name = "singleR.main")
rewildingproject <- AddMetaData(object = rewildingproject, metadata= predfine.RPmatrix$labels, col.name = "singleR.fine")
View(rewildingproject@meta.data)

saveRDS(rewildingproject, file = "RewildingProject.rds")
rewildingproject <- readRDS(file = "RewildingProject.rds")

Idents(rewildingproject) <- "singleR.main"
rewildingproject <- RunUMAP(rewildingproject, dims = 1:20)
DimPlot(rewildingproject, label = TRUE, reduction = "umap", repel = TRUE)

Idents(rewildingproject) <- "singleR.fine"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
View(rewildingproject@meta.data)

features <- "T cells (T.8Nve)"

Idents(rewildingproject) <- "batch"
FeaturePlot(rewildingproject, features = "batch")

cells.highlight = 
DimPlot(rewildingproject, label = TRUE, reduction = "umap", repel = TRUE, cells.highlight = "Batch1")

library(tidyverse)

umap_batch = rewildingproject@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(batch = rewildingproject@meta.data$batch)

ggplot(umap_batch, aes(x=UMAP_1, y=UMAP_2, color=batch)) + geom_point() + 
  scale_color_manual(values=c("Batch1" = "blue", 
                              "Batch2" = "red"))

Idents(rewildingproject) <- 'singleR.fine'
umap_fine = rewildingproject@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(fine = rewildingproject@meta.data$singleR.fine)

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.8MEM.OT1.D45.LISOVA)" = "blue", 
                              "T cells (T.8Mem)" = "blue",
                              "T cells (T.8MEM)" = "blue",
                              "T cells (T.8MEM.OT1.D106.VSVOVA)" = "blue",
                              "T cells (T.4Mem)" = "red",
                              "T cells (T.4MEM)" = "red",
                              "T cells (T.4MEM44H62L)" = "red"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "blue", 
                              "T cells (T.CD4.24H)" = "blue",
                              "T cells (T.CD4.48H)" = "blue",
                              "T cells (T.CD4.5H)" = "blue",
                              "T cells (T.CD4.CTR)" = "blue",
                              "T cells (T.CD4+TESTDB)" = "blue",
                              "T cells (T.CD4+TESTNA)" = "blue",
                              "T cells (T.CD4CONTROL)" = "blue",
                              "T cells (T.CD4TESTCJ)" = "red",
                              "T cells (T.CD4TESTJS)" = "blue"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "blue", 
                              "T cells (T.CD4.24H)" = "blue",
                              "T cells (T.CD4.48H)" = "blue",
                              "T cells (T.CD4.5H)" = "blue",
                              "T cells (T.CD4.CTR)" = "blue",
                              "T cells (T.CD4+TESTDB)" = "blue",
                              "T cells (T.CD4+TESTNA)" = "blue",
                              "T cells (T.CD4CONTROL)" = "blue",
                              "T cells (T.CD4TESTCJ)" = "red",
                              "T cells (T.CD4TESTJS)" = "blue",
                              "T cells (T.CD8.1H)" = "green2",
                              "T cells (T.CD8.48H)" = "green2",
                              "T cells (T.CD8.5H)" = "green2",
                              "T cells (T.CD8.CTR)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD8.1H)" = "green2",
                              "T cells (T.CD8.48H)" = "green2",
                              "T cells (T.CD8.5H)" = "green2",
                              "T cells (T.CD8.CTR)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "blue", 
                              "T cells (T.CD4.24H)" = "blue",
                              "T cells (T.CD4.48H)" = "blue",
                              "T cells (T.CD4.5H)" = "blue",
                              "T cells (T.CD4.CTR)" = "blue",
                              "T cells (T.CD4+TESTDB)" = "blue",
                              "T cells (T.CD4+TESTNA)" = "blue",
                              "T cells (T.CD4CONTROL)" = "blue",
                              "T cells (T.CD4TESTCJ)" = "red",
                              "T cells (T.CD4TESTJS)" = "blue",
                              "T cells (T.CD8.1H)" = "green2",
                              "T cells (T.CD8.48H)" = "green2",
                              "T cells (T.CD8.5H)" = "green2",
                              "T cells (T.CD8.CTR)" = "green2",
                              "T cells (T.Tregs)" = "gold"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "red", 
                              "T cells (T.CD4.24H)" = "red",
                              "T cells (T.CD4.48H)" = "red",
                              "T cells (T.CD4.5H)" = "red",
                              "T cells (T.CD4.CTR)" = "red",
                              "T cells (T.CD4+TESTDB)" = "red",
                              "T cells (T.CD4+TESTNA)" = "red",
                              "T cells (T.CD4CONTROL)" = "red",
                              "T cells (T.CD4TESTCJ)" = "red",
                              "T cells (T.CD4TESTJS)" = "red",
                              "T cells (T.CD8.1H)" = "blue",
                              "T cells (T.CD8.48H)" = "blue",
                              "T cells (T.CD8.5H)" = "blue",
                              "T cells (T.CD8.CTR)" = "blue"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "red", 
                              "T cells (T.CD4.24H)" = "red",
                              "T cells (T.CD4.48H)" = "red",
                              "T cells (T.CD4.5H)" = "red",
                              "T cells (T.CD4.CTR)" = "red",
                              "T cells (T.CD4+TESTDB)" = "red",
                              "T cells (T.CD4+TESTNA)" = "red",
                              "T cells (T.CD4CONTROL)" = "red",
                              "T cells (T.CD4TESTCJ)" = "gold",
                              "T cells (T.CD4TESTJS)" = "red"))


ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.Tregs)" = "red",
                              "T cells (T.8EFF.OT1.12HR.LISOVA)" = "blue",
                              "T cells (T.8EFF.OT1.24HR.LISOVA)" = "blue",
                              "T cells (T.8EFF.OT1.48HR.LISOVA)" = "blue",
                              "T cells (T.8EFF.OT1.D10LIS)" = "blue",
                              "T cells (T.8EFF.OT1.D5.VSVOVA)" = "blue",
                              "T cells (T.8EFF.OT1LISO)" = "blue",
                              "T cells (T.8EFF.TBET-.OT1LISOVA)" = "blue"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("B cells (B.FrF)" = "red",
                              "B cells (B.GC)" = "blue",
                              "B cells (B.MEM)" = "gold",
                              "B cells (B.MZ)" = "green2"))
                              
ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("B cells (B.T2)" = "red",
                              "B cells (B.T3))" = "blue",
                              "B cells (B1a)" = "green2",
                              "B cells (B1A)" = "green2",
                              "B cells (B1b)" = "gold"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.4.PLN)" = "red", 
                              "T cells (T.4EFF49D+11A+.D8.LCMV)" = "gold",
                              "T cells (T.4FP3-)" = "blue",
                              "T cells (T.4FP3+25+)" = "blue"))


ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("DC (DC)" = "blue", 
                              "DC (DC.8-4-11B-)" = "blue",
                              "DC (DC.8-4-11B+)" = "blue",
                              "DC (DC.8+)" = "blue",
                              "Monocytes (MO.6C+II-)" = "red",
                              "Monocytes (MO.6C-II+)" = "red",
                              "Tgd (Tgd.VG2+)" = "gold",
                              "Tgd (Tgd.vg2-.act)" = "gold",
                              "Tgd (Tgd.vg2-)" = "gold",
                              "Tgd (Tgd.mat.VG1+VD6+)" = "gold",
                              "Tgd (Tgd.mat.VG2+)" = "gold",
                              "Tgd (Tgd.mat.vg3)" = "gold",
                              "Macrophages (MF.103-11B+.SALM3)" = "green2",
                              "Macrophages (MF.RP)" = "green2",
                              "Macrophages (MF)" = "green2",
                              "Macrophages (MF.II-480HI)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("Basophils (BA)" = "blue", 
                              "ILC (ILC3.LTI.4+)" = "red",
                              "ILC (ILC2)" = "red",
                              "ILC (ILC3.LTI.CD4-)" = "red",
                              "ILC (ILC3.LTI.CD4+)" = "red",
                              "Monocytes (MO.6+2+)" = "gold",
                              "Monocytes (MO.6C-II+)" = "gold",
                              "Monocytes (MO.6C+II-)" = "gold",
                              "Monocytes (MO.6C-IIINT)" = "gold",
                              "Neutrophils (GN.ARTH)" = "green2",
                              "Neutrophils (GN.Thio)" = "green2",
                              "Neutrophils (GN.URAC)" = "green2",
                              "Neutrophils (GN)" = "green2",
                              "Stem cells (GMP)" = "violet",
                              "Stem cells (MLP)" = "violet",
                              "Stem cells (proB.CLP)" = "violet",
                              "Stem cells (SC.CDP)" = "violet",
                              "Stem cells (SC.LT34F)" = "violet",
                              "Stem cells (SC.MDP)" = "violet",
                              "Stem cells (SC.MEP)" = "violet"))

#Rename cell clusters
View(rewildingproject@meta.data)
Idents(rewildingproject) <- "seurat_clusters"
new.cluster.ids <- c("CD4", "B.Fo", "CD8", "B.CD83", "NKT", "NK", "T.IFN", "Treg", 
                     "B.pro", "T.Tnfrsf9", "B.GC.Aicda", "CD8.Eff", "CD4.Nr4a1", "MO.MF", "B.unk", "DC.Fscn1", "B.Dctpp1", "DC.Plbd1",
                     "B.Mzb1.1", "B.IFN", "B.Mzb1.2", 
                    "T.lncRNA", "pDC")
names(new.cluster.ids) <- levels(rewildingproject)
rewildingproject <- RenameIdents(rewildingproject, new.cluster.ids)
DimPlot(rewildingproject, label = TRUE, reduction = "umap")
rewildingproject[["CellType"]] <- Idents(object = rewildingproject)
View(rewildingproject@meta\
     
     
     
     
     
     
     
     
     
     
     
     
     .data)




















Idents(rewildingproject) <- "CellType"
Idents(rewildingproject) <- "RNA"
CD4Tcells <- subset(rewildingproject, idents = c("CD4", "CD4.Nr4a1"))
FeaturePlot(CD4Tcells, features = "Il4", pt.size = 1.0, order = TRUE)
FeaturePlot(rewildingproject, features = "Il4", pt.size = 1.0, order = TRUE)
FeaturePlot(rewildingproject, features = "Il13", pt.size = 1, order = TRUE, split.by = "strain")
DotPlot(object = rewildingproject, features = "Il5", split.by = "Block")
Idents(rewildingproject) <- "environmentgroup"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "infectionstatus"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "strain"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "Block"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")
Idents(rewildingproject) <- "orig.ident"
DimPlot(rewildingproject, label = FALSE, reduction = "umap")




#scPCA
Idents(rewildingproject) <- "HTO_maxID"
PCA <- sc_pca(rewildingproject, "CellType", "hash.ID")
ggsave("PCA.pdf", PCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/Rewilding Project/Block 1 + Block 2/")

#Block PCA
blockPCA <- sc_pca(rewildingproject, "CellType", "hash.ID", "Block")
ggsave("BlockPCA.pdf", blockPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/Rewilding Project/Block 1 + Block 2/")

#Environment PCA
environmentPCA <- sc_pca(rewildingproject, "CellType", "hash.ID", "environmentgroup")
ggsave("environmentPCA.pdf", environmentPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/Rewilding Project/Block 1 + Block 2/")

#Infection PCA
infectionPCA <- sc_pca(rewildingproject, "CellType", "hash.ID", "infectionstatus")
ggsave("infectionPCA.pdf", infectionPCA[[1]], height=10, width=20,path = "/Users/howardng/Documents/Rewilding Project/Block 1 + Block 2/")

#Genotype PCA
genotypePCA <- sc_pca(rewildingproject, "CellType", "hash.ID", "strain")
ggsave("genotypePCA.pdf", genotypePCA[[1]], height=10, width=20, path = "/Users/howardng/Documents/Rewilding Project/Block 1 + Block 2/")

View(rewilding.combined@meta.data)
rewilding.combined[["CellType"]] <- Idents(object = rewilding.combined)

saveRDS(rewildingproject, file = "rewildingprojectobject.rds")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install dittoSeq
BiocManager::install("dittoSeq")
library(dittoSeq)

dittoBarPlot(
  object = rewildingproject,
  var = "CellType",
  group.by = "environmentgroup",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

dittoBarPlot(
  object = rewildingproject,
  var = "CellType",
  group.by = "strain",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

dittoBarPlot(
  object = rewildingproject,
  var = "CellType",
  group.by = "infectionstatus",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

dittoBarPlot(
  object = rewildingproject,
  var = "CellType",
  group.by = "Block",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

dittoBarPlot(
  object = rewildingproject,
  var = "CellType",
  group.by = "HTO_maxID",
  do.hover = T,
  theme =theme_classic( base_size = 8
  ))