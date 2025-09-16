library(Seurat)
library(clusterProfiler)

#Load data and create seurat object
Batch1.data <- Read10X(data.dir = "/Users/vandana/Desktop/Vandana/xiaohuaintestinalsample/Batch1/filtered_feature_bc_matrix")
Batch1 <- CreateSeuratObject(counts = Batch1.data, project = "Batch1", min.cells = 3, min.features = 200)
Batch1[["percent.mt"]] <- PercentageFeatureSet(Batch1, pattern = "^MT-")
VlnPlot(Batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Perform QC and filter dataset
Batch1 <- subset(Batch1, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & percent.mt < 25)
#Normalize data
Batch1 <- NormalizeData(Batch1)
Batch1 <- FindVariableFeatures(Batch1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Batch1)
#Scale data
Batch1 <- ScaleData(Batch1, features = all.genes)
#Run PCA, dimensionality reduction and clustering
Batch1 <- RunPCA(Batch1, features = VariableFeatures(object = Batch1))
ElbowPlot(Batch1)
Batch1 <- FindNeighbors(Batch1, dims = 1:15)
Batch1 <- FindClusters(Batch1, resolution = 0.5)
Batch1 <- RunUMAP(Batch1, dims = 1:15)
DimPlot(Batch1, reduction = "umap")

saveRDS(Batch1, "/athena/chenlab/scratch/jjv4001/Qianglab/Batch1.rds" )

#Load data and create seurat object
MIT.data <- Read10X(data.dir = "/Users/vandana/Desktop/Vandana/xiaohuaintestinalsample/MIT/filtered_feature_bc_matrix")
MIT <- CreateSeuratObject(counts = MIT.data, project = "Batch1", min.cells = 3, min.features = 200)
MIT[["percent.mt"]] <- PercentageFeatureSet(MIT, pattern = "^MT-")
VlnPlot(MIT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Perform QC and filter dataset
MIT <- subset(MIT, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & percent.mt < 25)
#Normalize data
MIT <- NormalizeData(MIT)
MIT <- FindVariableFeatures(MIT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MIT)
#Scale data
MIT <- ScaleData(MIT, features = all.genes)
#Run PCA, dimensionality reduction and clustering
MIT <- RunPCA(MIT, features = VariableFeatures(object = MIT))
ElbowPlot(MIT)
MIT <- FindNeighbors(MIT, dims = 1:15)
MIT <- FindClusters(MIT, resolution = 0.5)
MIT <- RunUMAP(MIT, dims = 1:15)
DimPlot(MIT, reduction = "umap")

saveRDS(MIT, "/athena/chenlab/scratch/jjv4001/Qianglab/MIT.rds" )

#Load data and create seurat object
Human.data <- Read10X(data.dir = "/Users/vandana/Desktop/Vandana/xiaohuaintestinalsample/Human/filtered_feature_bc_matrix")
Human <- CreateSeuratObject(counts = Human.data, project = "Batch1", min.cells = 3, min.features = 200)
Human[["percent.mt"]] <- PercentageFeatureSet(Human, pattern = "^MT-")
VlnPlot(Human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Perform QC and filter dataset
Human <- subset(Human, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 20)
#Normalize data
Human <- NormalizeData(Human)
Human <- FindVariableFeatures(Human, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Human)
#Scale data
Human <- ScaleData(Human, features = all.genes)
#Run PCA, dimensionality reduction and clustering
Human <- RunPCA(Human, features = VariableFeatures(object = Human))
ElbowPlot(Human)
Human <- FindNeighbors(Human, dims = 1:15)
Human <- FindClusters(Human, resolution = 0.5)
Human <- RunUMAP(Human, dims = 1:15)
DimPlot(Human, reduction = "umap")

saveRDS(Human, "/athena/chenlab/scratch/jjv4001/Qianglab/Human.rds" )


#Create function for assigning cell clusters based on marker gene expression
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


#Load data from iPSC organoids 
Batch1<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/Batch1.rds")
MIT<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/MIT.rds")
Batch1$SampleType<-rep("Batch1")
MIT$SampleType<-rep("MIT")
#merge data from iPSC organoids
organoids<-merge(Batch1, MIT)
organoids<- JoinLayers(organoids)
#normalize data from iPSC organoids
organoids <- NormalizeData(organoids)
organoids <- FindVariableFeatures(organoids)
#scale data from iPSC organoids
organoids <- ScaleData(organoids)
#Run PCA, dimensionality reduction and clustering
organoids <- RunPCA(organoids)
organoids <- RunUMAP(organoids, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
#integrate samples using harmony
organoids<-RunHarmony(object=organoids, group.by.vars="SampleType")
organoids <- RunUMAP(organoids, dims = 1:15, reduction = "harmony", reduction.name = "umap.integrated")

#Define markers for assignment
marker_lists <- list(
  Epithelial = c("EPCAM", "CDX2"),
  Mesenchymal = c("VIM","FBLN5")
)

#assign epithelial/mesenchymal identity clusters 
organoids<- assign_celltypes_by_score(organoids, marker_lists, min_score = 0.3)
DimPlot(organoids, group.by = "predicted_celltype", label = TRUE, repel = TRUE)
Idents(organoids)<-organoids$predicted_celltype

#subset epithelial cells from iPSC organoids
Epithelial<-subset(organoids, idents="Epithelial")
#normalize data
Epithelial <- NormalizeData(Epithelial)
Epithelial <- FindVariableFeatures(Epithelial)
#Scale data
Epithelial <- ScaleData(Epithelial)
#Run PCA, dimensionality reduction and clustering
Epithelial <- RunPCA(Epithelial)
Epithelial <- RunUMAP(Epithelial, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
#assign epithelial cell types 
marker_lists <- list(
  Stem = c("LGR5","RGMB","BMI1"),
  Enterocytes = c("CEACAM5","CEACAM6","FABP1"),
  TA = c("TOP2A","MKI67","PCNA"),
  Goblet = c("MUC2","MUC5AC","SPINK4"),
  EEC  = c("PAX6", "CHGB","CHGA")
  
)
Epithelial<- assign_celltypes_by_score(Epithelial, marker_lists, min_score = 0.3)
saveRDS(Epithelial, "/athena/chenlab/scratch/jjv4001/Qianglab/iPSCorganoidsassigned.rds")

#Load human organoid data
Human<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/Human.rds")
#Define markers for assignment
marker_lists <- list(
  Epithelial = c("EPCAM", "CDX2"),
  Mesenchymal = c("VIM","FBLN5")
)

#assign epithelial/mesenchymal identity clusters 
Human<- assign_celltypes_by_score(Human, marker_lists, min_score = 0.3)
DimPlot(organoids, group.by = "predicted_celltype", label = TRUE, repel = TRUE)
Idents(Human)<-organoids$predicted_celltype

#subset epithelial cells from iPSC organoids
Epithelial1<-subset(Human, idents="Epithelial")
#normalize data
Epithelial1 <- NormalizeData(Epithelial)
Epithelial1 <- FindVariableFeatures(Epithelial1)
#Scale data
Epithelial1 <- ScaleData(Epithelial1)
#Run PCA, dimensionality reduction and clustering
Epithelial1 <- RunPCA(Epithelial1)
Epithelial1 <- RunUMAP(Epithelial1, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
#assign epithelial cell types 
marker_lists <- list(
  Stem = c("LGR5","RGMB","BMI1"),
  Enterocytes = c("CEACAM5","CEACAM6","FABP1"),
  TA = c("TOP2A","MKI67","PCNA"),
  Goblet = c("MUC2","MUC5AC","SPINK4"),
  EEC  = c("PAX6", "CHGB","CHGA")
  
)
Epithelial1<- assign_celltypes_by_score(Epithelial1, marker_lists, min_score = 0.3)
saveRDS(Epithelial1, "/athena/chenlab/scratch/jjv4001/Qianglab/humanorganoidsassigned.rds")

#Load all annotated smaples
iPSC<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/iPSCorganoidsassigned.rds")
human<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/humanorganoidsassigned.rds")
intetine<-readRDS("/athena/chenlab/scratch/jjv4001/Qianglab/intestinewithoutpediatric.rds")
#specify identity and sampletype
iPSC$SampleType<-rep("iPSC")
iPSC$Identity<-Idents(iPSC)
human$SampleType<-rep("human")
human$Identity<-Idents(human)
intetine$Identity<-Idents(intetine)
#merge hPSC and primary organoids with human fetal and adult intestinal samples
all<-merge(x=intetine, y=c(iPSC, human), add.cell.ids = c("Tissue","iPSCorganoids","Humanorganoids"), project = "intestinalorganoids",
           merge.data = TRUE)
all<- JoinLayers(all)
#normalize merged data
all <- NormalizeData(all)
all <- FindVariableFeatures(all)
#scale merged data
all <- ScaleData(all)
all <- RunPCA(all)
all <- RunUMAP(all, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
#save merged data
saveRDS(all, file=paste("allsamplesmerged.rds", sep=""))
library(harmony)
#integrate merged data
all<-RunHarmony(object=all, group.by.vars="SampleType")
all <- RunUMAP(all, dims = 1:15, reduction = "harmony", reduction.name = "umap.integrated")
Idents(all)<-all$Identity
#save merged data
saveRDS(all, file=paste("allsamplesintegrated1.rds", sep=""))

#find top expressed markers in the adult
Idents(all)<-all$SampleType
adult.markers<-FindMarkers(all, ident.1 ="adult" , ident.2 = "fetal")
adult.markers.sig<-subset(adult.markers, p_val_adj<0.05)
#remove immunoglobulin genes 
deg_filtered <- adult.markers.sig[!grepl("^IGLV|^IGKV|^IGL|^IGK|^IGH", rownames(adult.markers.sig), ignore.case = TRUE), ]
#extract top 100 gene markers
top_100 <- deg_filtered[order(-deg_filtered$avg_log2FC), ][1:100, ]
#find average expression of all genes in each SampleType
avg_expr<-data.frame(AverageExpression(all))
write.csv(avg_expr, "/athena/chenlab/scratch/jjv4001/Qianglab/averageexpressionall.csv")

#create new csv files for average expression of all genes in iPSC and primary human organoids in one folder
#and run GSEA
inputgenes<-rownames(top_100)
#customGSEA, input genes as a string, e.g. c("INS", "PDX1", etc.)
customGSEA<-function(inputgenes, outputfile){
  genes<-inputgenes
  nlength<-length(genes)
  #TERM2GENE, prepare custom geneset using user input genes, geneset requires a 2 column data frame with pathway/term on the left and genes on the right 
  gene_sets<-data.frame(Term=rep("Custom Gene List", nlength), Gene=genes)
  #Load our DE lists for GSEA for each hormone into a single list of dataframes
  directory_path <- "/athena/chenlab/scratch/jjv4001/Qianglab/deseq"
  csv_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)
  list_of_data_frames <- lapply(csv_files, read.csv)
  names(list_of_data_frames) <- tools::file_path_sans_ext(basename(csv_files))
  #Convert the dataframes into named vectors of decreasing log2foldchange, with genes as names
  
  #convert_to_named_vector <- function(df, value_var, name_var) {
  #  vec <- setNames(df[[value_var]], df[[name_var]])
  #  sorted_vec <- sort(vec, decreasing = TRUE)
  #  return(sorted_vec)
  #}
  convert_to_named_vector <- function(df, value_var, name_var) {
    df <- df[!is.na(df[[value_var]]) & !is.na(df[[name_var]]), ]  # Remove NAs
    df <- df[!duplicated(df[[name_var]]), ]  # Remove duplicates
    vec <- setNames(df[[value_var]], df[[name_var]])
    sorted_vec <- sort(vec, decreasing = TRUE)
    return(vec)
  }
  named_vectors_list <- lapply(list_of_data_frames, convert_to_named_vector, value_var = "avg_expr", name_var = "gene_symbol")
  order_named_vector <- function(named_vector) {
    named_vector[order(-named_vector)]
  }
  ordered_named_vectors_list <- lapply(named_vectors_list, order_named_vector)
  #Run GSEA for each hormone against custom gene list, custom gene list is supplied in TERM2GENE and genelist supplied consists of ordered named numeric vectors and create a list of the GSEA results
  run_gsea <- function(gene_list, gene_sets) {
    gsea_result <- GSEA(geneList = gene_list, TERM2GENE = gene_sets, minGSSize=5, maxGSSize=500, pvalueCutoff = 1, pAdjustMethod="BH", verbose=TRUE, nPerm=1000)
    return(gsea_result)
  }
  
  
  gsea_results_list <- lapply(ordered_named_vectors_list, run_gsea, gene_sets = gene_sets)
  for (i in seq_along(gsea_results_list)) {
    cat("GSEA Results for Vector", names(ordered_named_vectors_list)[i], ":\n")
    print(gsea_results_list[[i]])
  }
  
  filtered_gsea_results_list <- Filter(function(x) {
    dims <- dim(x)
    !(dims[1] == 0)
  }, gsea_results_list)
  #convert GSEA list to a dataframe containing GSEA results from all hormones 
  gsea_to_df <- function(gsea_result, identifier) {
    df <- as.data.frame(gsea_result)
    df$perturbation_id <- identifier  # Add an identifier column to track the source of each result
    return(df)
  }
  gsea_dfs <- lapply(seq_along(filtered_gsea_results_list), function(i) {
    gsea_to_df(filtered_gsea_results_list[[i]], names(filtered_gsea_results_list)[i])
  })
  combined_gsea_df <- do.call(rbind, gsea_dfs)
  write.csv(combined_gsea_df,outputfile)}

customGSEA(inputgenes, "GSEA.csv")