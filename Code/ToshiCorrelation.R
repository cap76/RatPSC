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

#Load mouse
#mouse_data0 <- readRDS(file=paste(saveext,"/MouseMerge",".rds",sep=""))
#mouse_data <- subset(mouse_data,idents=c("E7.5_Gut","E7.5_NotLabelledscent_mesoderm","E7.5_Notochord","E7.5_Pharyngeal_mesoderm"),invert=TRUE)
#mouse_data0 <- subset(mouse_data,idents=c("E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E4.5_NotLabelled"),invert=TRUE)

Ratmeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQCRAT.txt",sep="\t",header = T)
Mousemeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQMOUSE.txt",sep="\t",header = T)
key_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_sample_metadata.txt",sep="\t",header = T, row.names=1)

raw_counts1<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_Rat.txt",sep="\t",header = T, row.names=1)
raw_counts2<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_mousePE.txt",sep="\t",header = T, row.names=1)
counts_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_rna_counts_gn.txt", sep="\t", header=T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)

rownames(raw_counts1) <- toupper(rownames(raw_counts1))
rownames(raw_counts2) <- toupper(rownames(raw_counts2))
rownames(counts_mouse2) <- toupper(rownames(counts_mouse2))
rownames(raw_counts3) <- toupper(rownames(raw_counts3))

common_genes <- intersect(intersect(intersect(rownames(raw_counts1),rownames(raw_counts2)),rownames(counts_mouse2)),rownames(raw_counts3))

#Plot and visualise mouse

rat_data1 <- CreateSeuratObject(counts = raw_counts1[common_genes,], assay = "RNA",min.cells = 0, min.features = 0)
Idents(rat_data1) <- Ratmeta$ID
rat_data1$Batch <- "Rat"
rat_data1 <- NormalizeData(rat_data1, verbose = FALSE)
rat_data1 <- FindVariableFeatures(rat_data1, selection.method = "vst", nfeatures = 3000)

mouse_data1 <- CreateSeuratObject(counts = raw_counts2[common_genes,which(Mousemeta$QC>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data1) <- Mousemeta$ID[which(Mousemeta$QC>0)]
mouse_data1$Batch <- "Mouse1"
mouse_data1 <- NormalizeData(mouse_data1, verbose = FALSE)
mouse_data1 <- FindVariableFeatures(mouse_data1, selection.method = "vst", nfeatures = 3000)

mouse_data2 <- CreateSeuratObject(counts = counts_mouse2[common_genes,which(key_mouse2$pass_rnaQC==TRUE)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data2) <- key_mouse2$lineage10x_3[which(key_mouse2$pass_rnaQC==TRUE)]
mouse_data2 <- subset(mouse_data2, subset = nFeature_RNA > 0)
mouse_data2 <- NormalizeData(mouse_data2, verbose = FALSE)
mouse_data2 <- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 3000)


mergeddata <- merge(rat_data1, y = c(mouse_data1,mouse_data2), project = "merged")


avexp  <- AverageExpression(object = mergeddata) #, use.scale = TRUE)

l1 <- c("N11_Epi","N12_Epi","Rat_ES", "rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","rESC_PVC_HT","rESC_PVC_WT")
l2 <- c("FS1","FS2","AF","AFX")
l3 <- c("E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E7.5_Epiblast","E6.5_Primitive_Streak")

C1 <- cor(  log2(avexp$RNA[,l1] +1),log2(avexp$RNA[,l2] +1) )
C2 <- cor(  log2(avexp$RNA[,l1] +1),log2(avexp$RNA[,l3] +1) )

mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1_allgenes.pdf",sep=""))


mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2_allgenes.pdf",sep=""))



#####Now do more


#Tediously add all the folders can shortcut this later when we have finalised everything.

BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D_ALLKEY.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)
rownames(raw_counts3) <- toupper(rownames(raw_counts3))
mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)

isinvit <- BS$IsEmb
labs <- BS$ID2
labs2 <- BS$PosID
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
raw_counts3<- raw_counts3[,which(isinvit>0)]
mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
mouse_data$species <- "mouse"
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)
Idents(mouse_data) <- labs
mouse_data$ID <- labs
mouse_data$ID2 <- labs2


Idents(mouse_data) <- mouse_data$ID2

avexp2  <- AverageExpression(object = mouse_data) #, use.scale = TRUE)

C3 <- cor(  log2(avexp$RNA[common_genes,l1] +1),log2(avexp2$RNA[common_genes,] +1) )


write.csv(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs.csv",sep=""))
write.csv(as.data.frame(Idents(mouse_data)), file=paste(saveext,"/Header.csv",sep=""))

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



