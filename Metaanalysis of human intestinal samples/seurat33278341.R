#Load Seurat
library(Seurat)
# Specify markers for visualization
genes <- c("EPCAM", "CDX2","COL1A1","COL6A1","TUBB2B", 
           "ELAVL4","HBA1","HBB","PECAM1","GSN",
           "CCL3","CD74")
# PMID of the dataset
id <- "PMID33278341"
message(paste("sample:", id))

#Load the respective matrix files for run and merge the matrices
dirs <- dir(path="./", pattern="run_count")
seurat <- lapply(dirs, function(dir){
  message(dir)
  dir <- paste("./", dir, "/outs/filtered_feature_bc_matrix/", sep="")
  seurat <- Read10X(dir)
  seurat <- CreateSeuratObject(counts=seurat, project=gsub("/outs/filtered_feature_bc_matrix/", "", gsub("./run_count_", "", dir)))
  seurat})
seurat <- merge(seurat[[1]], y=seurat[2:length(seurat)], project=id, add.cell.ids=gsub("run_count_", "", dirs))
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern="^MT-")

# data preparation, filter cells with nFeature_RNA<100 and nFeature_RNA > 10000, percent.mt > 40)
seurat <- subset(seurat, subset=nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 40)

# normalize data
seurat <- NormalizeData(seurat, normalization.method="LogNormalize", scale.factor=10000)
# identification of highly variable features 
seurat <- FindVariableFeatures(seurat, selection.method="vst", nfeatures=400)
# scale data
seurat <- ScaleData(seurat, features=rownames(seurat))
# run linear dimensional reduction
seurat <- RunPCA(seurat , features=VariableFeatures(object=seurat))
# run non-linear dimensional reduction
seurat <- RunUMAP(seurat, dims=1:10, n.components=3)
seurat@misc$umap3d <- seurat@reductions$umap
seurat <- RunUMAP(seurat, dims=1:10, n.components=2)
seurat@misc$umap2d <- seurat@reductions$umap
# cluster cells
seurat <- FindNeighbors(seurat, dims=1:10)
pdf(paste(id, "_umap_clusters_iter_1.pdf", sep=""), width=5, height=5)
for (resolution in seq(0, 1, 0.05)){
  g <- FindClusters(seurat, resolution=resolution)
  g <- DimPlot(g, reduction="umap", label=TRUE, raster=TRUE)
  print(g)}
dev.off()

# visualization
pdf(paste(id, "_umap_QC_iter_1.pdf", sep=""), width=9, height=8)
FeaturePlot(seurat, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=TRUE)
dev.off()

pdf(paste(id, "_umap_markers_iter_1.pdf", sep=""), width=16, height=48)
FeaturePlot(seurat, features=genes, raster=TRUE)
dev.off()

pdf(paste(id, "_umap_group_iter_1.pdf", sep=""))
DimPlot(seurat, reduction="umap", group.by="orig.ident", label=TRUE, raster=TRUE)
dev.off()

# save Seurat analysis
saveRDS(seurat, file=paste(id, "_iter_1_seurat.rds", sep=""))


####################################################
##########          iteration #1          ##########
####################################################

# Load saved Seurat analysis
merged.gut<-readRDS("/home/fs01/jjv4001/PMID33278341_iter_1_seurat.rds")

# cluster cells
merged.gut <- FindClusters(merged.gut, resolution=0.15)
# assign cluster identities
cells <- structure(c("Mesenchymal", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Mesenchymal","Immune cells","Immune cells","Epithelial","Neuronal","Endothelial","Epithelial","RBCs"), names=levels(merged.gut))
merged.gut <- RenameIdents(merged.gut, cells)
merged.gut$Cell_Type <- as.character(merged.gut@active.ident)
#umap visualization
pdf(paste("PMID33278341_umap_cluster_iter_6.pdf", sep=""))
DimPlot(merged.gut, reduction="umap", label=TRUE)
dev.off()
# subset Epithelial cells
merged.gut<-subset(merged.gut, idents="Epithelial")
# save Seurat analysis
saveRDS(merged.gut, file=paste("PMID33278341_iter_2_seurat.rds", sep=""))
