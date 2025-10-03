library(Seurat)
#Load xenium data for all tissue samples
Section2<-LoadXenium(
    "/Volumes/Jeya/output-XETG00085__0057724__17241_2__20250512__155203",
     fov = "fov",
     assay = "Xenium",
     mols.qv.threshold = 20,
     cell.centroids = TRUE,
     molecule.coordinates = TRUE,
     segmentations = NULL,
     flip.xy = FALSE
   )

Section6<-LoadXenium(
     "/Volumes/Jeya/output-XETG00085__0057712__17241_6__20250512__155203",
     fov = "fov",
     assay = "Xenium",
     mols.qv.threshold = 20,
     cell.centroids = TRUE,
     molecule.coordinates = TRUE,
     segmentations = NULL,
     flip.xy = FALSE
   )

Section5<-LoadXenium(
    "/Volumes/Jeya/output-XETG00085__0057724__17241_5__20250512__155203",
    fov = "fov",
     assay = "Xenium",
    mols.qv.threshold = 20,
    cell.centroids = TRUE,
     molecule.coordinates = TRUE,
    segmentations = NULL,
     flip.xy = FALSE
   )

Section7<-LoadXenium(
     "/Volumes/Jeya/output-XETG00085__0057712__17241_7__20250512__155203",
     fov = "fov",
     assay = "Xenium",
     mols.qv.threshold = 20,
     cell.centroids = TRUE,
     molecule.coordinates = TRUE,
     segmentations = NULL,
     flip.xy = FALSE
  )

Section1<-LoadXenium(
     "/Volumes/Jeya/output-XETG00085__0057724__17241_1__20250512__155203",
     fov = "fov",
     assay = "Xenium",
     mols.qv.threshold = 20,
     cell.centroids = TRUE,
     molecule.coordinates = TRUE,
     segmentations = NULL,
     flip.xy = FALSE
   )

#merge sections from liver samples and save liver sample data
Section1$sample<-rep("Section1")
Section2$sample<-rep("Section2")
liver<-merge(x=Section1, y=Section2, cell.ids=c("Section1","Section2"))
liver <- subset(liver, subset = nCount_Xenium > 0)
saveRDS(liver, "/Users/vandana/Desktop/Vandana/Qianglab/liver.rds")

#save intestinal sample data
Section7 <- subset(Section7, subset = nCount_Xenium > 0)
saveRDS(Section7, "/Users/vandana/Desktop/Vandana/Qianglab/intestine.rds")

#merge sections from brain samples and save brain sample data
Section5$sample<-rep("Section5")
Section6$sample<-rep("Section6")
brain<-merge(x=Section5, y=Section6, cell.ids=c("Section5","Section6"))
brain <- subset(brain, subset = nCount_Xenium > 0)
saveRDS(brain, "/Users/vandana/Desktop/Vandana/Qianglab/brain.rds")

#Identify human gene panel with variance less than or equal to 0.2 between different tissue types
brain<-JoinLayers(brain)
brainexpr <- GetAssayData(brain, assay = "Xenium", slot = "counts")
liver<-JoinLayers(liver)
liverexpr <- GetAssayData(liver, assay = "Xenium", slot = "counts")
intestineexpr <- GetAssayData(intestine, assay = "Xenium", slot = "counts")
humangenes<-read.csv("/Users/vandana/Desktop/Vandana/Qianglab/Humangenes.csv")
human_genes<-humangenes$genes
brain_matrix_subset <- brainexpr[rownames(brainexpr) %in% human_genes, ]
liver_matrix_subset <- liverexpr[rownames(liverexpr) %in% human_genes, ]
intestine_matrix_subset <- intestineexpr[rownames(intestineexpr) %in% human_genes, ]
avg_expr <- data.frame(
  brain = rowMeans(brain_matrix_subset),
  liver = rowMeans(liver_matrix_subset),
  intestine = rowMeans(intestine_matrix_subset)
)
cv <- apply(avg_expr, 1, function(x) sd(x) / mean(x))
uniform_genes <- rownames(avg_expr)[cv < 0.2] 

#Define human gene panel with variance less than or equal to 0.2 between different tissue types
human_genes<-c("APOBEC3B","CARD8","DEFB103B","DIRAS3","GYG2","INSL4","LPA","NLRP11","SERPINA3","SHOX","SLC2A4RG","TNFRSF10C","TSPY1","ZNF469","MG833229.1-Yeast")

#Load save xenium data for each tissue type 
liver<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/liver.rds")
intestine<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/intestine.rds")
brain<-readRDS("/Users/vandana/Desktop/Vandana/Qianglab/brain.rds")
  
# Get expression matrix for human markers for liver sample, intestinal sample, and brain samples
human_counts <- GetAssayData(liver, assay = "Xenium", slot = "counts")[human_genes, ]
human_counts1 <- GetAssayData(intestine, assay = "Xenium", slot = "counts")[human_genes, ]
human_counts2 <- GetAssayData(brain, assay = "Xenium", slot = "counts")[human_genes, ]

# Identify cells expressing ≥1 human genes from human gene panel identified above for liver, intestinal and brain samples
cells_with_human <- colnames(human_counts)[colSums(human_counts > 0) >= 1]
cells_with_human1 <- colnames(human_counts1)[colSums(human_counts > 0) >= 1]
cells_with_human2 <- colnames(human_counts2)[colSums(human_counts > 0) >= 1]

# Subset Seurat object with human cells for liver sample, intestinal sample and brain samples
liverHuman <- subset(liver, cells = cells_with_human)
saveRDS(liverHuman, "/Users/vandana/Desktop/Vandana/Qianglab/liverhumangenes.rds")

intestineHuman <- subset(intestine, cells = cells_with_human1)
saveRDS(intestineHuman, "/Users/vandana/Desktop/Vandana/Qianglab/intestinehumangenes.rds")

brainHuman <- subset(brain, cells = cells_with_human2)
saveRDS(brainHuman, "/Users/vandana/Desktop/Vandana/Qianglab/brainhumangenes.rds")

#Normalize intestinal section data, perform dimensional reduction and assign cell identities 
intestineHuman <- SCTransform(intestineHuman, assay = "Xenium")
intestineHuman <- RunPCA(intestineHuman, assay = "SCT", verbose = FALSE)
intestineHuman <- FindNeighbors(intestineHuman, reduction = "pca", dims = 1:30)
intestineHuman <- FindClusters(intestineHuman, resolution = 10.0)
intestineHuman <- RunUMAP(intestineHuman, reduction = "pca", dims = 1:30)

new.cluster.ids <- c("Enterocytes","Stem cells","Stem cells","Enterocytes","Goblet cells","Enterocytes",
                     "TA cells","Enterocytes","Stem cells","Enterocytes","Stem cells",
                     "TA cells","EECs","Goblet cells","TA cells","Goblet cells",
                     "Enterocytes","Enterocytes","Enterocytes","Enterocytes","Enterocytes",
                     "Enterocytes","Goblet cells","Enterocytes","Enterocytes","Enterocytes",
                     "Enterocytes","Stem cells","Enterocytes","Enterocytes","Goblet cells",
                     "EECs","Stem cells","Enterocytes","Enterocytes","Goblet cells",
                     "EECs","Goblet cells","Goblet cells","Goblet cells","EECs",
                     "TA cells","Enterocytes","Goblet cells","Enterocytes","TA cells",
                     "Enterocytes","Stem cells","Stem cells","TA cells","Goblet cells",
                     "Enterocytes","Enterocytes")
names(new.cluster.ids) <- levels(intestineHuman)
intestineHuman <- RenameIdents(intestineHuman, new.cluster.ids)
#saveData
saveRDS(intestineHuman,"/Users/vandana/Desktop/Vandana/Qianglab/intestineassigned.rds")
            
#Normalize brain section data, perform dimensional reduction and assign cell identities 
brainHuman <- SCTransform(brainHuman, assay = "Xenium")
brainHuman <- RunPCA(brainHuman, assay = "SCT", verbose = FALSE)
brainHuman<-RunHarmony(object=brainHuman, group.by.vars="sample")
brainHuman <- RunUMAP(brainHuman, dims = 1:30, reduction="harmony")
brainHuman <- FindNeighbors(brainHuman, reduction = "pca", dims = 1:30)
brainHuman <- FindClusters(brainHuman, resolution = 10.0)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes","Glutamatergic Neurons",
                     "6","Glutamatergic Neurons","GABAergic Neurons","Glutamatergic Neurons","GABAergic Neurons",
                     "11","12","13","14","GABAergic Neurons",
                     "GABAergic Neurons","17","18","Glutamatergic Neurons","GABAergic Neurons",
                     "21","22","Retinal progenitor cells","GABAergic Neurons","25",
                     "26","27","Glutamatergic Neurons","29","30",
                     "31","32","33","34","35",
                     "36","37","38","39","40",
                     "41","42","43","44","45",
                     "46","47","48","49","GABAergic Neurons",
                     "51","52","53","54","Retinal progenitor cells",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "11","Radial glial cells","13","14",
                     "17","18",
                     "21","22","Retinal progenitor cells","25",
                     "26","27","29","Radial glial cells",
                     "31","32","Glutamatergic Neurons","Radial glial cells","35",
                     "36","37","38","39","40",
                     "41","42","43","44","45",
                     "46","47","48","49",
                     "51","52","53","54",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "11","Radial glial cells","13","14",
                     "17","18",
                     "21","22","Retinal progenitor cells","25",
                     "26","27","29",
                     "31","32","35",
                     "36","Radial glial cells","38","39","40",
                     "41","42","43","44","45",
                     "46","47","Radial glial cells","49",
                     "51","52","53","54",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "11","Radial glial cells","Glutamatergic Neurons","14",
                     "17","18",
                     "21","Glutamatergic Neurons","Retinal progenitor cells","25",
                     "26","27","29",
                     "31","32","35",
                     "36","38","39","40",
                     "41","42","43","44","45",
                     "46","47","Glutamatergic Neurons",
                     "51","52","53","54",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Glutamatergic Neurons","Radial glial cells","Glutamatergic Neurons",
                     "Glutamatergic Neurons","18",
                     "21","Retinal progenitor cells","25",
                     "26","27","29",
                     "31","32","35",
                     "36","38","39","40",
                     "41","42","43","44","45",
                     "46","47",
                     "51","52","53","54",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells","Glutamatergic Neurons",
                     "Glutamatergic Neurons","Retinal progenitor cells","25",
                     "26","Glutamatergic Neurons","29",
                     "31","32","35",
                     "36","38","39","40",
                     "41","42","43","44","45",
                     "46","47",
                     "51","52","53","54",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","25",
                     "26","29",
                     "31","32","35",
                     "36","Glutamatergic Neurons","39","40",
                     "41","42","Glutamatergic Neurons","44","45",
                     "46","47",
                     "51","52","53","Glutamatergic Neurons",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","25",
                     "26","29",
                     "31","32","35",
                     "36","39","40",
                     "41","42","44","GABAergic Neurons",
                     "46","GABAergic Neurons",
                     "51","52","53",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","25",
                     "26","29",
                     "31","32","35",
                     "36","GABAergic Neurons","GABAergic Neurons",
                     "41","42","44",
                     "46",
                     "51","GABAergic Neurons","53",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","25",
                     "26","29",
                     "31","32","35",
                     "36",
                     "Glutamatergic Neurons","42","44",
                     "GABAergic Neurons",
                     "51","53",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("0","1","2","Glutamatergic Neurons","Astrocytes",
                     "6","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","25",
                     "26","29",
                     "31","32","35",
                     "36",
                     "Astrocytes","44",
                     "51","Astrocytes",
                     "56")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)


new.cluster.ids <- c("GABAergic Neurons","GABAergic Neurons","2","Glutamatergic Neurons","Astrocytes",
                     "Immune cells","GABAergic Neurons",
                     "Radial glial cells",
                     "Retinal progenitor cells","Immune cells",
                     "Immune cells","Mesenchymal cells",
                     "Mesenchymal cells","Intermediate progenitors","Endothelial cells",
                     "Intermediate progenitors","Endothelial cells",
                     "Intermediate progenitors","Intermediate progenitors")

names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("GABAergic Neurons","Intermediate progenitors","Glutamatergic Neurons","Astrocytes",
                     "Immune cells","Radial glial cells",
                     "Retinal progenitor cells","Mesenchymal cells","Intermediate progenitors",
                     "Endothelial cells")
names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)


new.cluster.ids <- c("GABAergic Neurons","Intermediate progenitors","Glutamatergic Neurons","Astrocytes+Radial glial cells",
                     "Immune cells","GABAergic Neurons",
                     "Retinal progenitor cells","Mesenchymal cells",
                     "Endothelial cells")
names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("GABAergic Neurons","Pre-myelinating oligodendrocytes","Glutamatergic Neurons","Astrocytes+Radial glial cells",
                     "Immune cells",
                     "Retinal progenitor cells","Mesenchymal cells",
                     "Endothelial cells")
names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("GABAergic Neurons","Pre-myelinating oligodendrocytes","Glutamatergic Neurons","Astrocytes+Radial glial cells",
                     "Endothelial cells",
                     "Intermediate progenitors","infiltrating T cells",
                     "Microglia")
names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)

new.cluster.ids <- c("GABAergic Neurons","Pre-myelinating oligodendrocytes","Glutamatergic Neurons","Astrocytes+Radial glial cells",
                     "Endothelial cells",
                     "Intermediate progenitors","Non-neural cells",
                     "Non-neural cells")
names(new.cluster.ids) <- levels(brainHuman)
brainHuman <- RenameIdents(brainHuman, new.cluster.ids)
#saveData            
saveRDS(brainHuman,"/Users/vandana/Desktop/Vandana/Qianglab/brainhumanassigned.rds")
            
#Normalize liver section data, perform dimensional reduction and assign cell identities 
liverHuman <- SCTransform(liverHuman, assay = "Xenium")
liverHuman <- RunPCA(liverHuman, assay = "SCT", verbose = FALSE)
liverHuman<-RunHarmony(object=liverHuman, group.by.vars="samples")
liverHuman <- RunUMAP(liverHuman, dims = 1:30, reduction="harmony")
liverHuman <- FindNeighbors(liverHuman, reduction = "pca", dims = 1:30)
liverHuman <- FindClusters(liverHuman, resolution = 2.8)

new.cluster.ids <- c("Hepatocytes","Hepatocytes","Hepatocytes","Hepatocytes","Hepatocytes","Endothelial cells",
                     "Hepatocytes","Kupffer cells","Hepatocytes","Other cells","Hepatocytes",
                     "Hepatic stellate cells","Hepatocytes","Hepatocytes","Hepatocytes","Hepatocytes",
                     "Hepatocytes","Hepatic stellate cells","Hepatocytes","Hepatocytes","Kupffer cells",
                     "Hepatocytes","Hepatocytes","Kupffer cells","Hepatoblasts","Hepatocytes",
                     "Hepatocytes","Hepatocytes","Hepatoblasts","Hepatocytes","Hepatic stellate cells",
                     "Hepatic stellate cells","Hepatoblasts","Hepatoblasts")
names(new.cluster.ids) <- levels(liverHuman)
liverHuman <- RenameIdents(liverHuman, new.cluster.ids)

new.cluster.ids <- c("Hepatocytes","Endothelial cells",
                     "Hepatocytes","Other cells",
                     "Hepatic stellate cells","Hepatoblasts")
names(new.cluster.ids) <- levels(liverHuman)
liverHuman <- RenameIdents(liverHuman, new.cluster.ids)
#saveData            
saveRDS(liverHuman,"/Users/vandana/Desktop/Vandana/Qianglab/liverassigned.rds")
