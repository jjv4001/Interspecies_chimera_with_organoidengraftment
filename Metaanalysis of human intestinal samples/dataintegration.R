#Load Seurat
library(Seurat)
#Load seurat datasets corresponding to each study
PMID33290721<-readRDS("/home/fs01/jjv4001/PMID33290721_iter_2_seurat.rds")
PMID33406409<-readRDS("/home/fs01/jjv4001/PMID33406409_iter_2_seurat.rds")
PMID31753849<-readRDS("/home/fs01/jjv4001/PMID31753849_iter_1_seurat.rds")
PMID33278341<-readRDS("/home/fs01/jjv4001/PMID33378341_iter_2_seurat.rds")
PMID34497389<-readRDS("/home/fs01/jjv4001/PMID34497389_iter_2_seurat.rds")
PMID35176508<-readRDS("/home/fs01/jjv4001/PMID35176508_iter_1_seurat.rds")
#Assign dataset to each study
PMID33290721$dataset<-"PMID33290721"
PMID33406409$dataset<-"PMID33406409"
PMID31753849$dataset<-"PMID31753849"
PMID33278341$dataset<-"PMID33278341"
PMID34497389$dataset<-"PMID34497389"
PMID35176508$dataset<-"PMID35176508"

#Merge all seurat objects
merged.gut <- merge(PMID33290721, y = c(PMID33406409, PMID31753849,PMID33278341,PMID34497389,PMID35176508), add.cell.ids = c("PMID33290721", "PMID33406409","PMID31753849", "PMID33278341","PMID34497389","PMID35176508"))
#Join Layers
merged.gut <- JoinLayers(merged.gut)
#Save merged seurat object
saveRDS(merged.gut, file=paste("merged_iter_1_seurat.rds", sep=""))
#Split merged object by each dataset
merged.gut[["RNA"]] <- split(merged.gut[["RNA"]], f = merged.gut$dataset)
#Normalize dataset
merged.gut <- NormalizeData(merged.gut)
# identification of highly variable features 
merged.gut <- FindVariableFeatures(merged.gut)
# scale data
merged.gut <- ScaleData(merged.gut)
# run linear dimensional reduction
merged.gut <- RunPCA(merged.gut)
# cluster cells
merged.gut<- FindNeighbors(merged.gut, dims = 1:10, reduction = "pca")
merged.gut <- FindClusters(merged.gut, resolution = 2, cluster.name = "unintegrated_clusters")
# run non-linear dimensional reduction
merged.gut <- RunUMAP(merged.gut, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")
# visualization
pdf(paste("unintegrated_umap_clusters_iter_1.pdf", sep=""), width=5, height=5)
DimPlot(merged.gut, reduction = "umap.unintegrated", group.by = c("dataset"), label=TRUE)+NoLegend()
dev.off()

options(future.globals.maxSize = 3e+09)
#Integrate the datasets
merged.gut <- IntegrateLayers(object = merged.gut, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",verbose = FALSE)
#JoinLayers
merged.gut[["RNA"]] <- JoinLayers(merged.gut[["RNA"]])
# cluster cells
merged.gut <- FindNeighbors(merged.gut, reduction = "integrated.cca", dims = 1:10)
merged.gut <- FindClusters(merged.gut, resolution = 1)
# run non-linear dimensional reduction
merged.gut <- RunUMAP(merged.gut, dims = 1:10, reduction = "integrated.cca", reduction.name = "umap.integrated")
# visualization
pdf(paste("integrated_umap_clusters_iter_4.pdf", sep=""), width=5, height=5)
DimPlot(merged.gut, reduction = "umap.integrated", split.by = c("dataset"))
dev.off()
pdf(paste("integrated_umap_clusters_iter_5.pdf", sep=""), width=5, height=5)
for (resolution in seq(0, 1, 0.05)){
  g <- FindClusters(merged.gut, resolution=resolution)
  g <- DimPlot(g, reduction="umap.integrated", label=TRUE, raster=TRUE)
  print(g)}
dev.off()

# cluster cells
merged.gut <- FindClusters(merged.gut, resolution=0.5)
# assign cluster identities
cells <- structure(c("Enterocytes", "Enterocytes", "Enterocytes", "Enterocytes", "Enterocytes", "Enterocytes","Enterocytes","Enterocytes","Enterocytes","Stem Cells","TA","Enterocytes","Enterocytes","Goblet Cells","Enterocytes","Stem Cells","Stem Cells","Enterocytes","TA","TA","Enterocytes","Enterocytes","Enterocytes","Enterocytes","Enterocytes","EECs","Enterocytes","Goblet Cells","Enterocytes","Enterocytes","Enterocytes","Enterocytes"), names=levels(merged.gut))
merged.gut <- RenameIdents(merged.gut, cells)
merged.gut$Cell_Type <- as.character(merged.gut@active.ident)
# visualization
pdf(paste("integrated_umap_cluster_iter_6.pdf", sep=""))
DimPlot(merged.gut, reduction="umap", label=TRUE)
dev.off()
# save Seurat analysis
saveRDS(merged.gut, file=paste("integrated_iter_2_seurat.rds", sep=""))

library(readxl)

#Load metadata information regarding sample type
metadata<-read_xlsx("/athena/chenlab/scratch/jjv4001/Qianglab/metadata.xlsx")
merged.gut$SampleType<-metadata$SampleType
#omit pediatric samples from analysis
mergedgutwoutpediatric<-subset(merged.gut, SampleType="pediatric", invert=TRUE)
saveRDS(mergedgutwoutpediatric,"/athena/chenlab/scratch/jjv4001/Qianglab/intestinewithoutpediatric.rds")


####################################################
##########          iteration #3          ##########
####################################################
