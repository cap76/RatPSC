#load various packages: Seurat for doing integration, ggplot for plots, pheatmap for heatmaps
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")

#Make sure the seed is set for repeatabilty
set.seed(1)

#Create folder for saving output.
saveext = "Heatmaps"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))


#List of cells in the final dataset
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

#List of colours for the cells in the variable "CellType"
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


#Load rate metadata
Ratmeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQCRAT.txt",sep="\t",header = T)

#Load the mouse in vitro meta data
Mousemeta<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQMOUSE.txt",sep="\t",header = T)
Mousemeta2<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/allQCB2.txt",sep="\t",header = T)

#load the mouse in vivo referece
key_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_sample_metadata.txt",sep="\t",header = T, row.names=1)

#Load the rat count data
raw_counts1<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_Rat.txt",sep="\t",header = T, row.names=1)

#Load the moouse in vitro count data
raw_counts2<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_mousePE.txt",sep="\t",header = T, row.names=1)
raw_counts2b<-read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My\ Drive/Toshi/Data/featurecountsAllgene_mousePE_b2.txt",sep="\t",header = T, row.names=1)

#Load the mouse in vivo count data
counts_mouse2 <- read.table("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Mouse/MouseSC1/GSE121650_rna_counts_gn.txt", sep="\t", header=T, row.names=1)

#Load the spatial count data
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)

#Extract out genes in each dataset and find the common genes
rownames(raw_counts1) <- toupper(rownames(raw_counts1))
rownames(raw_counts2) <- toupper(rownames(raw_counts2))
rownames(raw_counts2b) <- toupper(rownames(raw_counts2b))
rownames(counts_mouse2) <- toupper(rownames(counts_mouse2))
rownames(raw_counts3) <- toupper(rownames(raw_counts3))
common_genes <- intersect(intersect(intersect(intersect(rownames(raw_counts1),rownames(raw_counts2)),rownames(counts_mouse2)),rownames(raw_counts3)),rownames(raw_counts2b))

#Create seurat object of the rat data
rat_data1 <- CreateSeuratObject(counts = raw_counts1[common_genes,], assay = "RNA",min.cells = 0, min.features = 0)
Idents(rat_data1) <- Ratmeta$ID #assigns cell type to the cells 
rat_data1$Batch <- "Rat" #Assigns a designation to the cells in this dataset
rat_data1 <- NormalizeData(rat_data1, verbose = FALSE) #Log CP10k normalisation
rat_data1 <- FindVariableFeatures(rat_data1, selection.method = "vst", nfeatures = 3000)

#Create seurat object of in vitro data
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

#Create seurate object of in vivo
mouse_data2 <- CreateSeuratObject(counts = counts_mouse2[common_genes,which(key_mouse2$pass_rnaQC==TRUE)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(mouse_data2) <- key_mouse2$lineage10x_3[which(key_mouse2$pass_rnaQC==TRUE)]
mouse_data2 <- subset(mouse_data2, subset = nFeature_RNA > 0)
mouse_data2 <- NormalizeData(mouse_data2, verbose = FALSE)
mouse_data2 <- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 3000)

#Remove none embryonic cells from the in vivo reference
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)

#Merge the in vitro datasets
mergeddata <- merge(mouse_data1, y = c(mouse_data1b,rat_data1), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"

#Now align in vivo with in vitro usig 2000 integration genes
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"

#Normalise (zero mean, unit variance) for PCA
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

#Do PCA
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
#Do UMAP
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#Get nearest neighbours
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#Get the correct colours for each cells. This is the tedious way to do it, better ways to do so
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(CellType==levels(Idents(mammal.combined))[i])
}
coluse <- CellCol[colind]

#Now we will give the cells a label (0,1,2) depending on what experimment they are in. This will be used to shape the points in the plot. 0 = in vivo, 1 = mouse, 2 = rate
cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)%in%c("ES","FS1","FS2","AXR1","AXR2","AF","AFX") )] <- 1
cluster_letters[which(Idents(mammal.combined)%in%c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi","N12_Epi"))] <- 2
cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined, metadata = cluster_letters,col.name = 'cell.orig')

#Plot the UMAP
DimPlot(mammal.combined,  cols = coluse,  shape.by = "cell.orig", reduction = "umap",   label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_UMAP_split7",".pdf",sep=""),width = 15, height = 8,limitsize = FALSE)
#Plot PCA
DimPlot(mammal.combined, cols = coluse, shape.by = "cell.orig", reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/All_PCA_split7",".pdf",sep=""),width = 15, height = 8,limitsize = FALSE)
#Get list of 2000 genes use for the integration
ListInt1 <- rownames(GetAssayData(mammal.combined))

#Save the Seurat object (we can use it to do the DM reduction later)
saveRDS(mammal.combined,file=paste("ToshiInt.rds"))

#Get average expression
avexpma  <- AverageExpression(object = mammal.combined) #, use.scale = TRUE)

l1 <- c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi")
l2 <- c("ES","FS1","FS2","AXR2","AXR1","AF","AFX")
l3 <- c("E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E7.5_Epiblast","E6.5_Primitive_Streak")

#Get the correlation between rat and mouse in vitro
C1 <- cor(  log2(avexpma$integrated[ListInt1,l1] +1),log2(avexpma$integrated[ListInt1,l2] +1) )

#Get the correlation between rat and in vivo
C2 <- cor(  log2(avexpma$integrated[ListInt1,l1] +1),log2(avexpma$integrated[ListInt1,l3] +1) )

#Plot the heatmaps
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1Int_Int2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3),  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2Int_Int2kgenes.pdf",sep=""))


C1 <- cor(  log2(avexpma$RNA[ListInt1,l1] +1),log2(avexpma$RNA[ListInt1,l2] +1) )

#Get the correlation between rat and in vivo
C2 <- cor(  log2(avexpma$RNA[ListInt1,l1] +1),log2(avexpma$RNA[ListInt1,l3] +1) )

#Plot the heatmaps
mat_breaks <- seq(0.65, .8, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C1,color =  redblue1(20), display_numbers = round(C1, digits = 3), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C1RNA_Int2kgenes.pdf",sep=""))
mat_breaks <- seq(0.75, .85, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
pheatmap(C2,color =  redblue1(20), display_numbers = round(C2, digits = 3),  border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/LineageHeatmap_C2RNA_Int2kgenes.pdf",sep=""))




#Now do spatial analysis

#Reread the spatial data (metadata and count)
BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D_ALLKEY.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Mouse3D/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)
rownames(raw_counts3) <- toupper(rownames(raw_counts3))

#Create seurat object for spatial data
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

#Get list of cells that correpoond to endoderm
subsID <- colnames(mouse_data)[-grep("VE",mouse_data$ID3)]
mouse_data <- subset(mouse_data, cells=subsID )
subsID <- colnames(mouse_data)[-grep("EA",mouse_data$ID2)]
mouse_data <- subset(mouse_data, cells=subsID )
subsID <- colnames(mouse_data)[-grep("EP",mouse_data$ID2)]
mouse_data <- subset(mouse_data, cells=subsID )

#Get the mouse in vivo reference and remove endoderm
mouse_data2_sub <- subset(mouse_data2,idents=c("E4.5_Primitive_endoderm","E7.5_Visceral_endoderm","E5.5_Visceral_endoderm","E6.5_Visceral_endoderm","E7.5_Gut","E7.5_Pharyngeal_mesoderm","E6.5_NotLabelled","E6.5_NotLabelledscent_mesoderm","E7.5_Notochord","E4.5_NotLabelled","E7.5_NotLabelledscent_mesoderm"),invert=TRUE)
mergeddata <- merge(mouse_data1, y = c(mouse_data1b,rat_data1), project = "merged")
mouse_data2_sub$species <- "In vivo"
mergeddata$species <- "In vitro"
mouse_data$species <- "Spatial"
#Idents(mouse_data) <- "Spatial"


#Integrate the spatial data, in vivo reference, and the in vitro models
mammal.anchors <- FindIntegrationAnchors(object.list = list(mergeddata, mouse_data2_sub,mouse_data), dims = 1:20, k.anchor = 10, anchor.features = 2000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

#Plot UMAP and PCA
DimPlot(mammal.combined,  reduction = "umap",   label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/All_UMAP_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE) 
DimPlot(mammal.combined, reduction = "pca",  label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/All_PCA_split6",".pdf",sep=""),width = 10, height = 8,limitsize = FALSE) 
ListInt1 <- rownames(GetAssayData(mammal.combined))


#Redo the above
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

l1 <- c("Rat_ES","rESC_PVC_HT","rESC_PVC_WT","rEpiSC_NTAV","rEpiSC_WT_21","rEpiSC_WT_23","N11_Epi")
l2 <- c("ES","FS1","FS2","AXR2","AXR1","AF","AFX")
l3 <- as.character(unique(Idents(mouse_data)))

#Average expression
avexp2  <- AverageExpression(object = mammal.combined) #, use.scale = TRUE)

#Get correlations betwee in vitroo cells and each spatial location and save as csv file. Note we need to change this to a tab delimited text fille for use in Matlab scripts later.
C3 <- cor(  log2(avexp2$integrated[,l1] +1),log2(avexp2$integrated[,l3] +1) )
write.csv(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs1.csv",sep=""))
write.table(as.data.frame(C3), file=paste(saveext,"/AverageExpSubs1.txt",sep=""), sep="\t")

C4 <- cor(  log2(avexp2$integrated[,l1] +1),log2(avexp2$integrated[,l2] +1) )
write.csv(as.data.frame(C4), file=paste(saveext,"/AverageExpSubs2.csv",sep=""))
write.table(as.data.frame(C4), file=paste(saveext,"/AverageExpSubs2.txt",sep=""), sep="\t")


#Now we have the basic plots and the AverageExpSubs1.csv and AverageExpSubs2.csv correlations we can run the Matlab script to align samples. To do so open these csv files and save as tab-delimited .txt files.
#Next open MATLAB and run "main_Toshi.m" (you may need to change some file paths)
