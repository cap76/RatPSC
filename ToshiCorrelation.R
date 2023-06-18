library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "Heatmaps"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))


CellType <- c("E4.5_Epiblast",
             "E5.0",
             "E5.5",
             "E5.5_Epiblast",
             "E6.0",
             "E6.5_Epiblast",
             "E6.5_ExE_ectoderm",
             "E6.5_Parietal_endoderm",
             "E6.5_Primitive_Streak",
             "E6.5_Rostral_neurectoderm",
             "E7.5_Anterior_Primitive_Streak",
             "E7.5_Caudal_epiblast",
             "E7.5_Embryonic_endoderm",
             "E7.5_Epiblast",
             "E7.5_ExE_mesoderm",
             "E7.5_Mesenchyme",
             "E7.5_Mixed_mesoderm",
             "E7.5_Paraxial_mesoderm",
             "E7.5_Primitive_Streak",
             "E7.5_Rostral_neurectoderm",
             "E7.5_Surface_ectoderm",
             "ES",
             "FS1",
             "FS2",
             "AXR1",
             "AXR2",
             "AF",
             "AFX",
             "Rat_ES",
             "rESC_PVC_HT",
             "rESC_PVC_WT",
             "rEpiSC_NTAV",
             "rEpiSC_WT_21",
             "rEpiSC_WT_23",
             "N11_Epi",
             "N12_Epi")

CellCol <- c("#2bb6b8",
            "#2bb6b8",
            "#0c9cf5",
            "#0c9cf5",
            "#0767da",
            "#0767da",
            "#5454a1",
            "#5a325c",
            "#242f60",
            "#737073",
            "#4c2c4d",
            "#022999",
            "#613862",
            "#0233bf",
            "#967700",
            "#574a16",
            "#7a6417",
            "#473d13",
            "#181a33",
            "#636163",
            "#1a0873",
            "#8856a7",
            "#4292c6",
            "#08519c",
            "#c994c7",
            "#df65b0",
            "#e7298a",
            "#ce1256",
            "#8856a7",
            "#c994c7",
            "#df65b0",
            "#e7298a",
            "#ce1256",
            "#980043",
            "#67001f",
            "#67001f")


#Load mouse
#mouse_data0 <- readRDS(file=paste(saveext,"/MouseMerge",".rds",sep=""))
#mouse_data <- subset(mouse_data,idents=c("E7.5_Gut","E7.5_NotLabelledscent_mesoderm","E7.5_Notochord","E7.5_Pharyngeal_mesoderm"),invert=TRUE)
#mouse_data0 <- subset(mouse_data,idents=c("E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E4.5_NotLabelled"),invert=TRUE)

Ratmeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQCRAT.txt",sep="\t",header = T)
Mousemeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQMOUSE.txt",sep="\t",header = T)
Mousemeta2<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQCB2.txt",sep="\t",header = T)
key_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_sample_metadata.txt",sep="\t",header = T, row.names=1)

raw_counts1<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_Rat.txt",sep="\t",header = T, row.names=1)
raw_counts2<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_mousePE.txt",sep="\t",header = T, row.names=1)
raw_counts2b<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_mousePE_b2.txt",sep="\t",header = T, row.names=1)
counts_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_rna_counts_gn.txt", sep="\t", header=T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)

rownames(raw_counts1) <- toupper(rownames(raw_counts1))
rownames(raw_counts2) <- toupper(rownames(raw_counts2))
rownames(raw_counts2b) <- toupper(rownames(raw_counts2b))
rownames(counts_mouse2) <- toupper(rownames(counts_mouse2))
rownames(raw_counts3) <- toupper(rownames(raw_counts3))

common_genes <- intersect(intersect(intersect(intersect(rownames(raw_counts1),rownames(raw_counts2)),rownames(counts_mouse2)),rownames(raw_counts3)),rownames(raw_counts2b))

#Plot and visualise mouse
rat_data1 <- CreateSeuratObject(counts = raw_counts1[common_genes,], assay = "RNA",min.cells = 0, min.features = 0)
Idents(rat_data1) <- Ratmeta$ID
rat_data1$Batch <- "Rat"
rat_data1 <- NormalizeData(rat_data1, verbose = FALSE)
rat_data1 <- FindVariableFeatures(rat_data1, selection.method = "vst", nfeatures = 3000)

#in vitro
mouse_data1 <- CreateSeuratObject(counts = raw_counts2[common_genes,which(Mousemeta$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data1) <- Mousemeta$ID[which(Mousemeta$QC>0)]
mouse_data1$Batch <- "Mouse1"
mouse_data1 <- NormalizeData(mouse_data1, verbose = FALSE)
mouse_data1 <- FindVariableFeatures(mouse_data1, selection.method = "vst", nfeatures = 3000)
mouse_data1b <- CreateSeuratObject(counts = raw_counts2b[common_genes,which(Mousemeta2$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data1b) <- Mousemeta2$ID[which(Mousemeta2$QC>0)]
mouse_data1b$Batch <- "Mouse2"
mouse_data1b <- NormalizeData(mouse_data1b, verbose = FALSE)
mouse_data1b <- FindVariableFeatures(mouse_data1b, selection.method = "vst", nfeatures = 3000)
#In vivo
mouse_data2 <- CreateSeuratObject(counts = counts_mouse2[common_genes,which(key_mouse2$pass_rnaQC==TRUE)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data2) <- key_mouse2$lineage10x_3[which(key_mouse2$pass_rnaQC==TRUE)]
mouse_data2 <- subset(mouse_data2, subset = nFeature_RNA > 0)
mouse_data2 <- NormalizeData(mouse_data2, verbose = FALSE)
mouse_data2 <- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 3000)

mergeddata <- mouse_data2 #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 2000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
List2K <- mergeddata2@assays$RNA@var.features

mergeddata <- mouse_data2 #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 4000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
List4K <- mergeddata2@assays$RNA@var.features

mergeddata <- mouse_data2 #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 6000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
List6K <- mergeddata2@assays$RNA@var.features

#Look at the emDisc
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata2_sub <- FindVariableFeatures(mouse_data2_sub, selection.method = "vst", nfeatures = 4000)
mergeddata2_sub <- ScaleData(mergeddata2_sub, verbose = FALSE)
mergeddata2_sub <- RunPCA(mergeddata2_sub, npcs = 20, verbose = FALSE)
mergeddata2_sub <- RunUMAP(mergeddata2_sub, reduction = "pca", dims = 1:20)
mergeddata2_sub <- FindNeighbors(mergeddata2_sub, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPSub_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCASub_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
List4KS <- mergeddata2_sub@assays$RNA@var.features

mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata2_sub <- FindVariableFeatures(mergeddata2_sub, selection.method = "vst", nfeatures = 2000)
mergeddata2_sub <- ScaleData(mergeddata2_sub, verbose = FALSE)
mergeddata2_sub <- RunPCA(mmergeddata2_sub, npcs = 20, verbose = FALSE)
mergeddata2_sub <- RunUMAP(mergeddata2_sub, reduction = "pca", dims = 1:20)
mergeddata2_sub <- FindNeighbors(mergeddata2_sub, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPSub_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCASub_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
List2KS <- mergeddata2_sub@assays$RNA@var.features

mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata2_sub <- FindVariableFeatures(mouse_data2_sub, selection.method = "vst", nfeatures = 6000)
mergeddata2_sub <- ScaleData(mergeddata2_sub, verbose = FALSE)
mergeddata2_sub <- RunPCA(mergeddata2_sub, npcs = 20, verbose = FALSE)
mergeddata2_sub <- RunUMAP(mergeddata2_sub, reduction = "pca", dims = 1:20)
mergeddata2_sub <- FindNeighbors(mergeddata2_sub, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPSub_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2_sub, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCASub_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
List6KS <- mergeddata2_sub@assays$RNA@var.features

#
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata2 <- merge(mouse_data2_sub, y = c(mouse_data1,mouse_data1b), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata2, selection.method = "vst", nfeatures = 4000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPMergSub_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCAMergSub_Type_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
List4KS <- mergeddata2@assays$RNA@var.features

mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata2 <- merge(mouse_data2_sub, y = c(mouse_data1,mouse_data1b), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata2, selection.method = "vst", nfeatures = 2000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPMergSub_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCAMergSub_Type_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
List2KS <- mergeddata2@assays$RNA@var.features


mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data2_sub, y = c(mouse_data1,mouse_data1b), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 6000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAPMergSub_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCAMergSub_Type_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
List6KS <- mergeddata2@assays$RNA@var.features



#
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data1, y = c(mouse_data1b), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined,  reduction = "umap",   label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_UMAP_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE)
DimPlot(mammal.combined, reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_PCA_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE)
ListInt1 <- rownames(GetAssayData(mammal.combined))

mouse_data2_sub <- subset(mouse_data2,idents=c("E7.5_Paraxial_mesoderm","E6.5_Parietal_endoderm","E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data1, y = c(mouse_data1b, rat_data1), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(CellType==levels(Idents(mammal.combined))[i])
}
coluse <- CellCol[colind]





cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)%in%c("ES","FS1","FS2","AXR1","AXR2","AF","AFX") )] <- 1
cluster_letters[which(Idents(mammal.combined)%in%c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi","N12_Epi"))] <- 2
cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined, metadata = cluster_letters,col.name = 'cell.orig')


DimPlot(mammal.combined,  cols = coluse,  shape.by = "cell.orig", reduction = "umap",   label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_UMAP_split7",".pdf",sep=""),width = 15, height = 8,limitsize = FALSE)
DimPlot(mammal.combined, cols = coluse, shape.by = "cell.orig", reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_PCA_split7",".pdf",sep=""),width = 15, height = 8,limitsize = FALSE)
ListInt1 <- rownames(GetAssayData(mammal.combined))

saveRDS(mammal.combined,file=paste("ToshiInt.rds"))

#

mergeddata <- merge(rat_data1, y = c(mouse_data1,mouse_data1b,mouse_data2), project = "merged")
avexp  <- AverageExpression(object = mergeddata) #, use.scale = TRUE)

l1 <- c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi","N12_Epi")
l2 <- c("ES","FS1","FS2","AXR2","AXR1","AF","AFX")
l3 <- c("E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E7.5_Epiblast","E6.5_Primitive_Streak")

C1 <- cor(  log2(avexp$RNA[,l1] +1),log2(avexp$RNA[,l2] +1) )
C2 <- cor(  log2(avexp$RNA[,l1] +1),log2(avexp$RNA[,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_allgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_allgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[List6K,l1] +1),log2(avexp$RNA[List6K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[List6K,l1] +1),log2(avexp$RNA[List6K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_6kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_6kgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[List4K,l1] +1),log2(avexp$RNA[List4K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[List4K,l1] +1),log2(avexp$RNA[List4K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_4kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_4kgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[List2K,l1] +1),log2(avexp$RNA[List2K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[List2K,l1] +1),log2(avexp$RNA[List2K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_2kgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[ListInt1,l1] +1),log2(avexp$RNA[ListInt1,l2] +1) )
C2 <- cor(  log2(avexp$RNA[ListInt1,l1] +1),log2(avexp$RNA[ListInt1,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3),  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_Int2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_Int2kgenes.pdf",sep=""))

avexpma  <- AverageExpression(object = mammal.combined) #, use.scale = TRUE)

l1 <- c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi","N12_Epi")
l2 <- c("ES","FS1","FS2","AXR2","AXR1","AF","AFX")
l3 <- c("E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E7.5_Epiblast","E6.5_Primitive_Streak")

C1 <- cor(  log2(avexpma$integrated[ListInt1,l1] +1),log2(avexpma$integrated[ListInt1,l2] +1) )
C2 <- cor(  log2(avexpma$integrated[ListInt1,l1] +1),log2(avexpma$integrated[ListInt1,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1Int_Int2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3),  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2Int_Int2kgenes.pdf",sep=""))


#####Now do more


#Tediously add all the folders can shortcut this later when we have finalised everything.

BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D_ALLKEY.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)
rownames(raw_counts3) <- toupper(rownames(raw_counts3))
mouse_data <- CreateSeuratObject(counts = raw_counts3[common_genes,], assay = "RNA",min.cells = 0, min.features = 0)

isinvit <- BS$IsEmb
labs <- BS$ID2
labs2 <- BS$PosID
labs3 <- BS$Desc
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
labs3 <- labs3[which(isinvit>0)]
raw_counts3<- raw_counts3[,which(isinvit>0)]
mouse_data <- CreateSeuratObject(counts = raw_counts3[common_genes,], assay = "RNA",min.cells = 0, min.features = 0)
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
mouse_data$species <- "mouse"
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)
Idents(mouse_data) <- labs
mouse_data$ID <- labs
mouse_data$ID2 <- labs2
mouse_data$ID3 <- labs3
Idents(mouse_data) <- mouse_data$ID2
subsID <- colnames(mouse_data)[-grep("VE",mouse_data$ID3)]
mouse_data <- subset(mouse_data, cells=subsID )
subsID <- colnames(mouse_data)[-grep("EA",mouse_data$ID2)]
mouse_data <- subset(mouse_data, cells=subsID )
subsID <- colnames(mouse_data)[-grep("EP",mouse_data$ID2)]
mouse_data <- subset(mouse_data, cells=subsID )





mergeddata <- mouse_data #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 2000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Spatial_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Spatial_Labs","_2k.pdf",sep=""),width = 20, height = 8,p)
ListSp2K <- mergeddata2@assays$RNA@var.features

mergeddata <- mouse_data #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 4000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Spatial_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Spatial_Labs","_4k.pdf",sep=""),width = 20, height = 8,p)
ListSp4K <- mergeddata2@assays$RNA@var.features

mergeddata <- mouse_data #merge(rat_data1, y = c(mouse_data2), project = "merged")
mergeddata2 <- FindVariableFeatures(mergeddata, selection.method = "vst", nfeatures = 6000)
mergeddata2 <- ScaleData(mergeddata2, verbose = FALSE)
mergeddata2 <- RunPCA(mergeddata2, npcs = 20, verbose = FALSE)
mergeddata2 <- RunUMAP(mergeddata2, reduction = "pca", dims = 1:20)
mergeddata2 <- FindNeighbors(mergeddata2, reduction = "pca", dims = 1:20)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Spatial_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
p <- DimPlot(mergeddata2, pt.size = 4, reduction = "pca", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Spatial_Labs","_6k.pdf",sep=""),width = 20, height = 8,p)
ListSp6K <- mergeddata2@assays$RNA@var.features


C1 <- cor(  log2(avexp$RNA[ListSp6K,l1] +1),log2(avexp$RNA[ListSp6K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[ListSp6K,l1] +1),log2(avexp$RNA[ListSp6K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_Sp6kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_Sp6kgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[ListSp4K,l1] +1),log2(avexp$RNA[ListSp4K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[ListSp4K,l1] +1),log2(avexp$RNA[ListSp4K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_Sp4kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_Sp4kgenes.pdf",sep=""))


C1 <- cor(  log2(avexp$RNA[ListSp2K,l1] +1),log2(avexp$RNA[ListSp2K,l2] +1) )
C2 <- cor(  log2(avexp$RNA[ListSp2K,l1] +1),log2(avexp$RNA[ListSp2K,l3] +1) )
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_Sp2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_Sp2kgenes.pdf",sep=""))



avexp2  <- AverageExpression(object = mouse_data) #, use.scale = TRUE)

C3 <- cor(  log2(avexp$RNA[common_genes,l1] +1),log2(avexp2$RNA[common_genes,] +1) )
write.csv(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs.csv",sep=""))
write.csv(as.data.frame(Idents(mouse_data)), file=paste(saveext,"/Header.csv",sep=""))

C4 <- cor(  log2(avexp$RNA[ListSp6K,l1] +1),log2(avexp2$RNA[ListSp6K,] +1) )
write.csv(as.data.frame(C4), file=paste(saveext,"/AverageExpSubs_6k.csv",sep=""))
write.csv(as.data.frame(Idents(mouse_data)), file=paste(saveext,"/Header.csv",sep=""))

C5 <- cor(  log2(avexp$RNA[ListSp4K,l1] +1),log2(avexp2$RNA[ListSp4K,] +1) )
write.csv(as.data.frame(C5), file=paste(saveext,"/AverageExpSubs_4k.csv",sep=""))


C6 <- cor(  log2(avexp$RNA[ListSp2K,l1] +1),log2(avexp2$RNA[ListSp2K,] +1) )
write.csv(as.data.frame(C6), file=paste(saveext,"/AverageExpSubs_2k.csv",sep=""))




#
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data1, y = c(mouse_data1b,rat_data1), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"
mouse_data$species <- "Spatial"
#Idents(mouse_data) <- "Spatial"
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub,mouse_data), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined,  reduction = "umap",   label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/All_UMAP_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE) 
DimPlot(mammal.combined, reduction = "pca",  label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/All_PCA_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE) 
ListInt1 <- rownames(GetAssayData(mammal.combined))

Idents(mouse_data) <- mouse_data$ID2
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data1, y = c(mouse_data1b,rat_data1), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"
mouse_data$species <- "Spatial"
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub,mouse_data), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

l1 <- c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi","N12_Epi")
l2 <- c("ES","FS1","FS2","AXR2","AXR1","AF","AFX")
l3 <- as.character(unique(Idents(mouse_data)))

avexp2  <- AverageExpression(object = mammal.combined) #, use.scale = TRUE)

C3 <- cor(  log2(avexp2$integrated[,l1] +1),log2(avexp2$integrated[,l3] +1) )
write.csv(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs1.txt",sep=""))


C4 <- cor(  log2(avexp2$integrated[,l1] +1),log2(avexp2$integrated[,l2] +1) )
write.csv(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs2.txt",sep=""))

#write.csv(as.data.frame(l3), file=paste(saveext,"/Header1.csv",sep=""))


#Idents(mouse_data) <- mouse_data$ID2
#AvExp <- AverageExpression(object = mouse_data)
#write.csv(as.data.frame(AvExp$RNA*100), file=paste(saveext,"/AverageExp.csv",sep=""))

#Dsubs <- AvExp$RNA[c("Sox2","Pou5f1","Nanog","Prdm14","Fbxo2","Tdgf1","Sox17","Mixl1","Hnf4a","T","Eomes","Mixl1","Foxa2","Otx2","Dppa2","Dppa4"),] #GetAssayData(object = AvExp , slot = 'data')[c("Sox2","Pou5f1","Nanog","Prdm14","Fbxo2","Tdgf1","Sox17","Mixl1","Hnf4a","T"),]
#write.csv(as.data.frame(Dsubs), file=paste(saveext,"/AverageExpSubs.csv",sep=""))

#mouse_data <- subset(mouse_data, idents = c("Myo_CS7","ReGland_CS5","ReGland_CS7","Gland_CS5","Gland_CS6","Gland_CS7","4-cell_CS2","8-cell_CS2","Zy_CS1","cMor_CS3"), invert = TRUE)
#mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
#mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 20000)
#mouse_data <- ScaleData(mouse_data, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
#mouse_data<- RunPCA(mouse_data, npcs = 30, verbose = FALSE)
#mouse_data <- RunUMAP(mouse_data, reduction = "pca", dims = 1:20)
#mouse_data <- RunTSNE(mouse_data, reduction = "pca", dims = 1:20)
#mouse_data <- FindNeighbors(mouse_data, reduction = "pca", dims = 1:20)
#mouse_data <- FindClusters(mouse_data, resolution = 0.5)
#mouse_data$Cl <- Idents(mouse_data)



