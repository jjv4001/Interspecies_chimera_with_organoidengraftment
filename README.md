# Interspecies_chimera_with_organoidengraftment

## Metaanalysis of human intestinal samples
seurat33209721.R, seurat33406409.R, seurat34497389.R, seurat33278341.R, seurat31753849.R, and seurat35176508.R provide codes for initial processing and Seurat analysis of scRNA data from human intestinal samples. Cell types from different germ layers are assigned and epithelial cells which are of interest are subset. dataintegration.R integrates epithelial cells from multiple scRNA datasets and assigns different epithelial cell type identities, removing pediatric samples from analysis. 

## GSEA
GSEA.R provides codes for initial processing and Seurat analysis of scRNA data from hiPSC organoids and primary human intestinal organoids. Cell types from different germ layers are assigned and epithelial cells which are of interest are subset, merged, and integrated with human fetal and adult intestinal tissue samples. Custom functions for cell type assignment and GSEA are also provided.

## Xenium
Xenium.R assigns different cell types in the chimeric brain, intestine, and liver using Seurat. 
