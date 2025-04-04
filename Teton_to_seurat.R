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
