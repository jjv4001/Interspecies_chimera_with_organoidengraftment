
#Load required packages
library(Seurat)
#Load matrix and create seurat object for intestine
intestine <- Read10X(data.dir = "D:/new data/xilis/NatBiotech revision/intestine/outs/raw_feature_bc_matrix")
intestine <- CreateSeuratObject(intestine, project="intestine")
#perform qc and filter out low quality cells
intestine[["percent.mt"]] <- PercentageFeatureSet(intestine, pattern = "^MT[-\\.]")
VlnPlot(intestine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(intestine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
library(patchwork)
plot1 <- FeatureScatter(intestine, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(intestine, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
intestine <- subset(intestine, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5)
#Normalize intestine data, find variable features, scale data
intestine <- NormalizeData(intestine)
intestine <- FindVariableFeatures(intestine, nfeatures = 5000)
intestine_top_features <- head(VariableFeatures(intestine), 20)
plot1 <- VariableFeaturePlot(intestine)
plot2 <- LabelPoints(plot = plot1, points = intestine_top_features, repel = TRUE)
plot1 + plot2
intestine <- ScaleData(intestine, vars.to.regress = c("nFeature_RNA", "percent.mt"))
#Run PCA
intestine <- RunPCA(intestine, npcs = 50)
ElbowPlot(intestine, ndims = ncol(Embeddings(intestine, "pca")))
PCHeatmap(intestine, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
#Run non-linear dimensional reduction
intestine <- RunUMAP(intestine, dims = 1:20)
#saveData
saveRDS(intestine, file="D:/new data/xilis/NatBiotech revision/intestine/intestine.rds")
#Load processed data
intestine<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/intestinexh.rds")
#Normalize intestine data, find variable features, scale data
intestine<-NormalizeData(intestine)
intestine<-FindVariableFeatures(intestine, selection.method = "vst", nfeatures=2000)
intestine<-ScaleData(intestine)
#Run PCA
intestine <- RunPCA(intestine, features = VariableFeatures(object = intestine))
#Run non-linear dimensional reduction
intestine <- FindNeighbors(intestine, dims = 1:10)
intestine <- FindClusters(intestine, resolution = 0.5)
intestine <- RunUMAP(intestine, dims = 1:10)
#Run function to assign cell types
assign_celltypes_by_score <- function(seurat_obj, marker_lists, score_name_prefix = "CellTypeScore", min_score = 0) {
  # 1. Run AddModuleScore on the gene sets
  seurat_obj <- AddModuleScore(seurat_obj, features = marker_lists, name = score_name_prefix)
  
  # 2. Extract score columns
  score_cols <- paste0(score_name_prefix, 1:length(marker_lists))
  scores <- seurat_obj@meta.data[, score_cols]
  
  # 3. Cell type names (in the same order as marker_lists)
  cell_types <- names(marker_lists)
  
  # 4. Assign each cell to cell type with max score
  max_idx <- max.col(scores, ties.method = "first")
  max_scores <- apply(scores, 1, max)
  
  predicted <- ifelse(max_scores >= min_score, cell_types[max_idx], "Unknown")
  
  # 5. Add predicted cell type to metadata
  seurat_obj$predicted_celltype <- predicted
  
  return(seurat_obj)
}

#Define different markers for the respective cell types to identify epithelial and stromal cells
marker_lists <- list(
  Epithelial = c("EPCAM", "CDX2"),
  Mesenchymal = c("VIM","FBLN5")
)

#Run the function on the seurat object and assign cell identities
intestine<- assign_celltypes_by_score(intestine, marker_lists, min_score = 0.3)
Idents(intestine)<-intestine$predicted_celltype

#filter out unknown cell types
intestine1<-subset(intestine, idents="Unknown", invert=TRUE)
#create identity column with different epithelial and stromal cell types
intestine1$Identity<-Idents(intestine1)
#Normalize subset intestine data, find variable features, scale data
intestine1<-NormalizeData(intestine1)
intestine1<-FindVariableFeatures(intestine1, selection.method = "vst", nfeatures=2000)
intestine1<-ScaleData(intestine1)
#Run PCA
intestine1 <- RunPCA(intestine1, features = VariableFeatures(object = intestine))
#Run non-linear dimensional reduction
intestine1 <- FindNeighbors(intestine1, dims = 1:10)
intestine1 <- FindClusters(intestine1, resolution = 0.5)
intestine1 <- RunUMAP(intestine1, dims = 1:10)
#Assign cell identities
Idents(intestine1)<-intestine1$Identity


#Define different markers for the respective epithelial cell types
marker_lists <- list(
  Stemcells = c("BMI1","MYC","LGR5"),
  TAcells = c("MKI67","EPHB2","TOP2A"),
  Enterocytes = c("ANPEP","SLC15A1"),
  Gobletcells = c("ABCA4","KCNMA1","MUC2","CLCA1"),
  EECs = c("FEV","TPH1")
)

#Run the function on the seurat object and assign cell identities
intestine<- assign_celltypes_by_score(intestine, marker_lists, min_score = 0.3)
Idents(intestine)<-intestine$predicted_celltype
#filter out unknown cell types
intestine1<-subset(intestine1, idents="Unknown", invert=TRUE)
#create identity column with different epithelial cell types
intestine1$Identity<-Idents(intestine1)
#Normalize subset intestine data, find variable features, scale data
intestine1<-NormalizeData(intestine1)
intestine1<-FindVariableFeatures(intestine1, selection.method = "vst", nfeatures=2000)
intestine1<-ScaleData(intestine1)
#Run PCA
intestine1 <- RunPCA(intestine1, features = VariableFeatures(object = intestine1))
#Run non-linear dimensional reduction
intestine1 <- FindNeighbors(intestine1, dims = 1:10)
intestine1 <- FindClusters(intestine1, resolution = 5.0)
intestine1 <- RunUMAP(intestine1, dims = 1:10)
#Assign cell identities
Idents(intestine1)<-intestine1$Identity
#saveData for assigned data
saveRDS(intestine1,"/Users/vandana/Desktop/Vandana/Qianglab/intestinejjv.rds")

#Load matrix and create seurat object for liver
liver <- Read10X(data.dir = "D:/new data/xilis/NatBiotech revision/liver/outs/raw_feature_bc_matrix")
liver <- CreateSeuratObject(liver, project="liver")
#perform qc and filter out low quality cells
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT[-\\.]")
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
library(patchwork)
plot1 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
liver <- subset(liver, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5)
#Normalize liver data, find variable features, scale data
liver <- NormalizeData(liver)
liver <- FindVariableFeatures(liver, nfeatures = 5000)
liver_top_features <- head(VariableFeatures(liver), 20)
plot1 <- VariableFeaturePlot(liver)
plot2 <- LabelPoints(plot = plot1, points = liver_top_features, repel = TRUE)
plot1 + plot2
liver <- ScaleData(liver, vars.to.regress = c("nFeature_RNA", "percent.mt"))
#Run PCA
liver <- RunPCA(liver, npcs = 50)
ElbowPlot(liver, ndims = ncol(Embeddings(liver, "pca")))
PCHeatmap(liver, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
#Run non-linear dimensional reduction
liver <- RunUMAP(liver, dims = 1:30)
#saveData
saveRDS(liver, file="D:/new data/xilis/NatBiotech revision/liver/liver.rds")
#Load processed data
liver<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/liverxh.rds")
#Normalize liver data, find variable features, scale data
liver<-NormalizeData(liver)
liver<-FindVariableFeatures(liver, selection.method = "vst", nfeatures=2000)
liver<-ScaleData(liver)
#Run PCA
liver <- RunPCA(liver, features = VariableFeatures(object = liver))
#Run non-linear dimensional reduction
liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.5)
liver <- RunUMAP(liver, dims = 1:10)

#Run function to assign cell types
assign_celltypes_by_score <- function(seurat_obj, marker_lists, score_name_prefix = "CellTypeScore", min_score = 0) {
  # 1. Run AddModuleScore on the gene sets
  seurat_obj <- AddModuleScore(seurat_obj, features = marker_lists, name = score_name_prefix)
  
  # 2. Extract score columns
  score_cols <- paste0(score_name_prefix, 1:length(marker_lists))
  scores <- seurat_obj@meta.data[, score_cols]
  
  # 3. Cell type names (in the same order as marker_lists)
  cell_types <- names(marker_lists)
  
  # 4. Assign each cell to cell type with max score
  max_idx <- max.col(scores, ties.method = "first")
  max_scores <- apply(scores, 1, max)
  
  predicted <- ifelse(max_scores >= min_score, cell_types[max_idx], "Unknown")
  
  # 5. Add predicted cell type to metadata
  seurat_obj$predicted_celltype <- predicted
  
  return(seurat_obj)
}

#Define different markers for the respective cell types
marker_lists <- list(
  Hepatoblasts = c("PROM1","HHEX"),
  Hepatocytes = c("PAH","ALB"),
  Hepaticstellatecells = c("TIMP1","RELN"),
  Kupffercells = c("TIMD4","CSF1R","CD5L"),
  Endothelialcells=c("KDR","FLT1")
)

#Run the function on the seurat object and assign cell identities
liver<- assign_celltypes_by_score(liver, marker_lists, min_score = 0.3)
Idents(liver)<-liver$predicted_celltype
#filter out unknown cell types
liver1<-subset(liver, idents="Unknown", invert=TRUE)
#create identity column with different liver cell types
liver1$Identity<-Idents(liver1)
#Normalize subset liver data, find variable features, scale data
liver1<-NormalizeData(liver1)
liver1<-FindVariableFeatures(liver1, selection.method = "vst", nfeatures=2000)
liver1<-ScaleData(liver1)
#Run PCA
liver1 <- RunPCA(liver1, features = VariableFeatures(object = liver1))
#Run non-linear dimensional reduction
liver1 <- FindNeighbors(liver1, dims = 1:10)
liver1 <- FindClusters(liver1, resolution = 5.0)
liver1 <- RunUMAP(liver1, dims = 1:10)
#Assign cell identities
Idents(liver1)<-liver1$Identity
#saveData for assigned data
saveRDS(liver1,"/Users/vandana/Desktop/Vandana/Qianglab/liverjjv.rds")

#Load matrix and create seurat object for brain
brain <- Read10X(data.dir = "D:/new data/xilis/NatBiotech revision/brain/outs/raw_feature_bc_matrix")
brain <- CreateSeuratObject(brain, project="brain")
#perform qc and filter out low quality cells
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT[-\\.]")
VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
library(patchwork)
plot1 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
brain <- subset(brain, subset = nFeature_RNA > 250 & nFeature_RNA < 8000 & percent.mt < 5)
#Normalize brain data, find variable features, scale data
brain<- NormalizeData(brain)
brain <- FindVariableFeatures(brain, nfeatures = 5000)
brain_top_features <- head(VariableFeatures(brain), 20)
plot1 <- VariableFeaturePlot(brain)
plot2 <- LabelPoints(plot = plot1, points = brain_top_features, repel = TRUE)
plot1 + plot2
brain <- ScaleData(brain, vars.to.regress = c("nFeature_RNA", "percent.mt"))
#Run PCA
brain <- RunPCA(brain, npcs = 50)
ElbowPlot(brain, ndims = ncol(Embeddings(brain, "pca")))
PCHeatmap(brain, dims = 1:30, cells = 500, balanced = TRUE, ncol = 4)
#Run non-linear dimensional reduction
brain <- RunUMAP(brain, dims = 1:30)
#saveData
saveRDS(brain, file="D:/new data/xilis/NatBiotech revision/brain/brain.rds")
#Load processed data
brain<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/brainxh.rds")
#Normalize brain data, find variable features, scale data
brain<-NormalizeData(brain)
brain<-FindVariableFeatures(brain, selection.method = "vst", nfeatures=2000)
brain<-ScaleData(brain)
#Run PCA
brain <- RunPCA(brain, features = VariableFeatures(object = brain))
#Run non-linear dimensional reduction
brain <- FindNeighbors(brain, dims = 1:10)
brain <- FindClusters(brain, resolution = 0.5)
brain <- RunUMAP(brain, dims = 1:10)
#Run function to assign cell types
assign_celltypes_by_score <- function(seurat_obj, marker_lists, score_name_prefix = "CellTypeScore", min_score = 0) {
  # 1. Run AddModuleScore on the gene sets
  seurat_obj <- AddModuleScore(seurat_obj, features = marker_lists, name = score_name_prefix)
  
  # 2. Extract score columns
  score_cols <- paste0(score_name_prefix, 1:length(marker_lists))
  scores <- seurat_obj@meta.data[, score_cols]
  
  # 3. Cell type names (in the same order as marker_lists)
  cell_types <- names(marker_lists)
  
  # 4. Assign each cell to cell type with max score
  max_idx <- max.col(scores, ties.method = "first")
  max_scores <- apply(scores, 1, max)
  
  predicted <- ifelse(max_scores >= min_score, cell_types[max_idx], "Unknown")
  
  # 5. Add predicted cell type to metadata
  seurat_obj$predicted_celltype <- predicted
  
  return(seurat_obj)
}
#Define different markers for the respective cell types
marker_lists <- list(
  Intermediateprogenitors = c("PAX6","NNAT","MYB"),
  GlutamatergicNeurons = c("NEUROD2","TBR1","GRIN2B","GLS"),
  GABAergicNeurons = c("GABBR1","GADD45B"),
  AstrocytesRadialglialcells = c("GLI3","SLC1A3","AQP4","SLC1A2"),
  Premyelinatingoligodendrocytes = c("CADM4","SOX10"),
  Endothelialcells=c("KDR","FLT1")
)
#Run the function on the seurat object and assign cell identities
brain<- assign_celltypes_by_score(brain, marker_lists, min_score = 0.3)
Idents(brain)<-brain$predicted_celltype
#filter out unknown cell types
brain1<-subset(brain, idents="Unknown", invert=TRUE)
#create identity column with different liver cell types
brain1$Identity<-Idents(brain1)
#Normalize subset brain data, find variable features, scale data
brain1<-NormalizeData(brain1)
brain1<-FindVariableFeatures(brain1, selection.method = "vst", nfeatures=2000)
brain1<-ScaleData(brain1)
#Run PCA
brain1 <- RunPCA(brain1, features = VariableFeatures(object = brain1))
#Run non-linear dimensional reduction
brain1 <- FindNeighbors(brain1, dims = 1:10)
brain1 <- FindClusters(brain1, resolution = 5.0)
brain1 <- RunUMAP(brain1, dims = 1:10)
#Assign cell identities
Idents(brain1)<-brain1$Identity
#saveData for assigned data
saveRDS(brain1,"/Users/vandana/Desktop/Vandana/Qianglab/brainjjv.rds")

