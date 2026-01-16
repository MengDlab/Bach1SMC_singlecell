library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(DoubletFinder)

# 1. Data Loading and Object Creation (Assuming two samples: WT and KO)
# Read WT data
data_wt <- Read10X(data.dir = "path/to/WT") 
seurat_wt <- CreateSeuratObject(counts = data_wt, project = "WT")
seurat_wt$group <- "WT" # Add group information

# Read KO data
data_ko <- Read10X(data.dir = "path/to/KO")
seurat_ko <- CreateSeuratObject(counts = data_ko, project = "KO")
seurat_ko$group <- "KO" # Add group information

# Store objects in a list for batch QC and doublet removal
sc_list <- list(WT = seurat_wt, KO = seurat_ko)

# ==============================================================================
for (i in names(sc_list)) {
  # --- 2.1 Calculate mitochondrial percentage ---
  sc_list[[i]][["percent.mt"]] <- PercentageFeatureSet(sc_list[[i]], pattern = "^MT-")
  
  # --- 2.2 Perform QC filtering ---
  # Description: gene expression < 200 or > 6000 or > 10% mito
  sc_list[[i]] <- subset(sc_list[[i]], 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 6000 & 
                           percent.mt < 10)
  
  sc_list[[i]] <- NormalizeData(sc_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  sc_list[[i]] <- FindVariableFeatures(sc_list[[i]], selection.method = "vst", nfeatures = 2000)
  sc_list[[i]] <- ScaleData(sc_list[[i]])
  sc_list[[i]] <- RunPCA(sc_list[[i]], verbose = FALSE)
  
  # --- Run DoubletFinder ---
  nExp_poi <- round(0.0075 * ncol(sc_list[[i]])) 
  
  # Execute DoubletFinder
  sc_list[[i]] <- doubletFinder_v3(sc_list[[i]], 
                                   PCs = 1:30, 
                                   pN = 0.25, 
                                   pK = 0.001, 
                                   nExp = nExp_poi, 
                                   reuse.pANN = FALSE, 
                                   sct = FALSE)
  
  df_col <- grep("DF.classifications", colnames(sc_list[[i]]@meta.data), value = TRUE)
  sc_list[[i]]@meta.data$Doublet_Class <- sc_list[[i]]@meta.data[[df_col]]
  
  # Filter to keep only Singlets
  sc_list[[i]] <- subset(sc_list[[i]], subset = Doublet_Class == "Singlet")
}

# Merge WT and KO data
sc_combined <- merge(sc_list[[1]], y = sc_list[[2]])
sc_combined <- NormalizeData(sc_combined, normalization.method = "LogNormalize", scale.factor = 10000)
sc_combined <- FindVariableFeatures(sc_combined, selection.method = "vst", nfeatures = 2000)
sc_combined <- ScaleData(sc_combined, verbose = FALSE)
sc_combined <- RunPCA(sc_combined, npcs = 50, verbose = FALSE)
sc_combined <- RunHarmony(sc_combined, group.by.vars = "group") 
sc_combined <- RunTSNE(sc_combined, 
                       reduction = "harmony", 
                       dims = 1:50)
sc_combined <- FindNeighbors(sc_combined, 
                             reduction = "harmony", 
                             dims = 1:50)
sc_combined <- FindClusters(sc_combined, resolution = 0.5)

# Plot t-SNE clustering
p1 <- DimPlot(sc_combined, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
# Plot t-SNE by group (WT vs KO)
p2 <- DimPlot(sc_combined, reduction = "tsne", group.by = "group") + ggtitle("Group")
print(p1 + p2)

# Save object
saveRDS(sc_combined, file = "seurat_processed_harmony.rds")

# Install DoRothEA and Viper
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("viper")
devtools::install_github("saezlab/dorothea")

# Install iTALK
devtools::install_github("Coolgenome/iTALK")

# ==============================================================================
# 1. Load packages and prepare data
# ==============================================================================
library(dorothea)
library(viper)
library(Seurat)
library(tidyverse)
library(pheatmap)

# sc_combined$cell_type <- Idents(sc_combined) # Ensure cell_type column exists

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs

# Filter regulons with high confidence
regulons <- regulons %>% 
  dplyr::filter(confidence %in% c("A", "B", "C"))

# Run Viper to infer TF activity
sc_combined <- run_viper(sc_combined, regulons,
                         options = list(method = "scale", minsize = 4, 
                                        eset.filter = FALSE, cores = 1, 
                                        verbose = FALSE))

DefaultAssay(sc_combined) <- "dorothea"
sc_combined <- ScaleData(sc_combined)

unique_cells <- "Smooth Muscle Cells" 
sub_obj <- subset(sc_combined, idents = unique_cells)

Idents(sub_obj) <- "group"

tf_markers <- FindMarkers(sub_obj, 
                          ident.1 = "KO", 
                          ident.2 = "WT", 
                          only.pos = FALSE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.1) # Threshold can be lower because activity scores have a different range than expression

print(head(tf_markers))
top_tfs <- rownames(tf_markers)[1:20]

# AverageExpression calculates the average value for each group
avg_activity <- AverageExpression(sub_obj, features = top_tfs, assays = "dorothea")
avg_activity <- avg_activity$dorothea
pheatmap(avg_activity, 
         scale = "row", 
         cluster_cols = FALSE, 
         main = paste0("TF Activity: WT vs KO in ", unique_cells))
DefaultAssay(sc_combined) <- "RNA"

# iTALK Analysis
library(iTALK)
library(dplyr)
prepare_italk_data <- function(seurat_obj, group_name) {
  sub_seurat <- subset(seurat_obj, subset = group == group_name)
  expr_mat <- as.data.frame(t(as.matrix(sub_seurat@assays$RNA@data)))
  expr_mat$cell_type <- Idents(sub_seurat)
  expr_mat$compare_group <- group_name
  
  return(expr_mat)
}

# Prepare data for WT and KO
italk_data_wt <- prepare_italk_data(sc_combined, "WT")
italk_data_ko <- prepare_italk_data(sc_combined, "KO")

comm_cats <- c('growth factor', 'cytokine', 'other')

# Analyze WT group
highly_exprs_genes_wt <- rawParse(italk_data_wt, top_genes = 50, stats = "mean")

comm_res_wt <- NULL
for(cat in comm_cats){
  res <- FindLR(highly_exprs_genes_wt, datatype='mean count', comm_type=cat)
  comm_res_wt <- rbind(comm_res_wt, res)
}

comm_res_wt <- comm_res_wt %>% arrange(desc(cell_from_mean_exprs * cell_to_mean_exprs))

# --- Analyze KO group ---
highly_exprs_genes_ko <- rawParse(italk_data_ko, top_genes = 50, stats = "mean")

comm_res_ko <- NULL
for(cat in comm_cats){
  res <- FindLR(highly_exprs_genes_ko, datatype='mean count', comm_type=cat)
  comm_res_ko <- rbind(comm_res_ko, res)
}
comm_res_ko <- comm_res_ko %>% arrange(desc(cell_from_mean_exprs * cell_to_mean_exprs))

par(mfrow=c(1,2)) # Set canvas to 1 row, 2 columns to display WT and KO simultaneously

# Plot WT
LRPlot(comm_res_wt[1:20, ], 
       datatype = 'mean count', 
       cell_col = structure(hue_pal()(length(unique(Idents(sc_combined)))), 
                            names = levels(Idents(sc_combined))), 
       link.arr.lwd = comm_res_wt$cell_from_mean_exprs[1:20], 
       link.arr.width = comm_res_wt$cell_to_mean_exprs[1:20])
title("WT Communication")

# Plot KO
LRPlot(comm_res_ko[1:20, ], 
       datatype = 'mean count', 
       cell_col = structure(hue_pal()(length(unique(Idents(sc_combined)))), 
                            names = levels(Idents(sc_combined))), 
       link.arr.lwd = comm_res_ko$cell_from_mean_exprs[1:20], 
       link.arr.width = comm_res_ko$cell_to_mean_exprs[1:20])
title("KO Communication")

