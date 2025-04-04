##load libraries
library(Seurat)
library(future)
plan("sequential")
#plan("multisession", workers = 10)
library(ggplot2)
library(cytoprofiling)

options(future.globals.maxSize = 1500 * 1024^2)


##create seurat object by cytoprofiling(I edited cytoprofiling_to_seurat() function. min.features and min.cells were changed to 0. This resolved the difference in the number of cells in RNA assay and Protein assay)
input_filename = "/your_path/RawCellStats.parquet"
aviti24_data <- load_cytoprofiling(input_filename)
#remove "GAPDH_1" and "NSB_" targets
aviti24_data <- aviti24_data[, !(grepl("^NSB", colnames(aviti24_data)) | grepl("^GAPDH_1", colnames(aviti24_data)))]
#RNA 
RNA_aviti24_data <- aviti24_data[, !endsWith(colnames(aviti24_data), ".B01")]
RNA_seurat_data <- cytoprofiling_to_seurat(RNA_aviti24_data)
#Protein
Protein_aviti24_data <- cbind(aviti24_data[,1:7], aviti24_data[, endsWith(colnames(aviti24_data), ".B01")])
Protein_seurat_data <- cytoprofiling_to_seurat(Protein_aviti24_data)
Protein_seurat_data <- RenameAssays(Protein_seurat_data, assay = "RNA", new.assay = "Protein")
#merge 2 modalities
RNA_seurat_data[["Protein"]] <- Protein_seurat_data[["Protein"]]
RNA_seurat_data@meta.data$nCount_Protein <- Protein_seurat_data@meta.data$nCount_Protein
RNA_seurat_data@meta.data$nFeature_Protein <- Protein_seurat_data@meta.data$nFeature_Protein
#save
saveRDS(RNA_seurat_data, file = "~/Teton_raw_counts.rds")

##after that, I filtered cells and normalized in Seurat. It means I didn't use filter_cells() and normalize_cytoprofiling() functions. (But I'm not sure this method is the best)
Teton_seurat_data <- RNA_seurat_data
rm(RNA_seurat_data)
Teton_seurat_data@meta.data$WellLabel <- factor(Teton_seurat_data@meta.data$WellLabel, 
                levels = c("HepG2_10k", "HepG2_8k", "PC9", "PC9_with_gef", "PC9_2", "PC9_3", "PC9_22", "PC9_23", "PC9_40", "PC9_41", "PC9-ER", "PC9-ER_with_gef"))
#filtering
VlnPlot(Teton_seurat_data, features = c("nCount_RNA", "nFeature_RNA", "nCount_Protein", "nFeature_Protein"), ncol = 2,
  log = TRUE, pt.size = 0, group.by = "WellLabel") + NoLegend()
Teton_seurat_data <- subset(
  x = Teton_seurat_data,
  subset = nCount_RNA > 20
) #> Removing 46280 cells missing data for vars requested

#normalization
#RNA
DefaultAssay(Teton_seurat_data) <- "RNA"
Teton_seurat_data <- NormalizeData(Teton_seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
Teton_seurat_data <- FindVariableFeatures(Teton_seurat_data, selection.method = "vst", nfeatures = 50)
all.genes <- rownames(Teton_seurat_data)
Teton_seurat_data <- ScaleData(Teton_seurat_data, features = all.genes)
Teton_seurat_data <- RunPCA(Teton_seurat_data, features = VariableFeatures(object = Teton_seurat_data), npcs=20)
ElbowPlot(Teton_seurat_data)#10 is good
Teton_seurat_data <- RunUMAP(Teton_seurat_data, dims = 1:10, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DimPlot(Teton_seurat_data, reduction = "umap.rna")

#Protein
DefaultAssay(Teton_seurat_data) <- "Protein"
Teton_seurat_data <- NormalizeData(object = Teton_seurat_data, normalization.method = "CLR", margin = 2)
all.genes <- rownames(Teton_seurat_data)
Teton_seurat_data <- ScaleData(Teton_seurat_data, features = all.genes)
VariableFeatures(Teton_seurat_data) <- rownames(Teton_seurat_data)  # since the panel is small, treat all features as variable.
Teton_seurat_data <- RunPCA(object = Teton_seurat_data, npcs = 20, verbose = FALSE, reduction.name = 'apca')
ElbowPlot(Teton_seurat_data)#10 looks good
Teton_seurat_data <- RunUMAP(object = Teton_seurat_data, reduction = 'apca', dims = 1:10, verbose = FALSE, reduction.name = 'umap.protein', reduction.key = 'proteinUMAP_')
DimPlot(Teton_seurat_data, reduction = "umap.protein")

#integrate by wnn
Teton_seurat_data <- FindMultiModalNeighbors(
  Teton_seurat_data, reduction.list = list("pca", "apca"), 
  dims.list = list(1:10, 1:10), modality.weight.name = "RNA.weight"
)
Teton_seurat_data <- RunUMAP(Teton_seurat_data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Teton_seurat_data <- FindClusters(Teton_seurat_data, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)#take quite long time
DimPlot(Teton_seurat_data, reduction = "wnn.umap", group.by = "WellLabel")
saveRDS(Teton_seurat_data, file = "~/Teton_finish_processed.rds")


