

### I will analyze the single cell RNA Seq data from the Zhao Q,.et al. (2022)
### article. I previously downloaded the relevant GSM files from GEO database.

### I set my current directory to the folder containing all the files:

setwd("D:/NEURON_RNA_SEQ_MAKALESİ/GSE192987_RAW")

data_dir <- getwd()

list.dirs(data_dir) ### these folders contain the relevant barcodes.tsv.gz, 
### features.tsv.gz, and matrix.mtx.gz files fro each samples obtained 
### by the application of cellRanger tool on the single cell RNA libraries
### prepared using 10x Genomics Chromium Single Cell 3′ Reagent 
### Kit v3 chemistry and sequenced on an Illumina NovaSeq6000 S4 sequencer.

### I will use the Seurat package for the analysis. After reading 10X chromium
### data and creating the Seurat objects, I will need to integrate them into
### a single Seurat object before moving to downstream analysis.

library(Seurat)

### we have each samples' cellRanger data inside specific folders. We will form
### specific direction paths for each sample folder 
### by using stri_join() function.
  
expression_matrix_1 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770202_vagus070519V1"))

expression_matrix_2 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770203_vagus070519V2"))
expression_matrix_3 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770204_vagus071819V1"))
expression_matrix_4 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770205_vagus071819V2"))
expression_matrix_5 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770206_vagus051818"))
expression_matrix_6 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770207_vagus010419"))
expression_matrix_7 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770208_vagusRC3"))
expression_matrix_8 <- Read10X(stringi::stri_join(data_dir,
                                                  "/GSM5770209_vagusRC4"))

### Form Seurat objects for each sample:

seurat_object_1 = CreateSeuratObject(counts = expression_matrix_1)
seurat_object_2 = CreateSeuratObject(counts = expression_matrix_2)
seurat_object_3 = CreateSeuratObject(counts = expression_matrix_3)
seurat_object_4 = CreateSeuratObject(counts = expression_matrix_4)
seurat_object_5 = CreateSeuratObject(counts = expression_matrix_5)
seurat_object_6 = CreateSeuratObject(counts = expression_matrix_6)
seurat_object_7 = CreateSeuratObject(counts = expression_matrix_7)
seurat_object_8 = CreateSeuratObject(counts = expression_matrix_8)

### Form a list:
 
seurat_object_lists <- list(seurat_object_1, seurat_object_2, seurat_object_3,
                            seurat_object_4, seurat_object_5, seurat_object_6, 
                            seurat_object_7, seurat_object_8)


#### normalize and identify variable features for each dataset independently:

my_list <- lapply(X = seurat_object_lists, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

###  we select the features that are repeatedly variable across datasets 
###  to be used for integration:

features <- SelectIntegrationFeatures(object.list = my_list)
 
   my_list <- lapply(X = my_list, FUN = function(x) {
       x <- ScaleData(x, features = features, verbose = FALSE)
       x <- RunPCA(x, features = features, verbose = FALSE)
    })

### We increase the memory limit value: 
   
memory.limit(size=56000) ### 

### We find the anchors for the integration of different Seurat objects.

### We will firstly use the fast integration method by using reciprocal PCA:

vagus.anchors <- FindIntegrationAnchors(object.list = my_list, 
                                        anchor.features = features, 
                                        reduction = "rpca")

### Another method for integration of Seurat objects:

### we identify the anchors using the FindIntegrationAnchors() function, 
### which takes a list of Seurat objects as input, and use these anchors 
### to integrate the two datasets together with IntegrateData(). 

vagus.anchors <- FindIntegrationAnchors(object.list = my_list, 
                                        anchor.features = features)

### After this step, following command was executed to form  
### an integrated object, 

### NOTE: The following command was run in a PC with 96 Gb Ram to avoid memory
### error.

vagus.combined <- IntegrateData(anchorset = vagus.anchors)

## we set the default assay in our integrated object:

DefaultAssay(vagus.combined) <- "integrated"


vagus.combined <- ScaleData(vagus.combined, verbose = FALSE)
vagus.combined <- RunPCA(vagus.combined, npcs = 30, verbose = FALSE)
vagus.combined <- RunUMAP(vagus.combined, reduction = "pca", dims = 1:30)
vagus.combined <- FindNeighbors(vagus.combined, reduction = "pca", dims = 1:30)
vagus.combined <- FindClusters(vagus.combined, resolution = 0.5)

p1 <- DimPlot(vagus.combined, reduction = "umap")

############################################################################

DefaultAssay(vagus.combined) <- "RNA"


DimPlot(vagus.combined, reduction = "umap")

### UMAP plot with labeled clusters:

DimPlot(vagus.combined, reduction = "umap", label = TRUE, repel = TRUE)

### Clusters that are rendered according to their expression of a specific gene:

FeaturePlot(vagus.combined, features = c("Slc17a6"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("Lypd6"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("P2rx2"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("Ptgdr"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("Mtpn"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("Prdm12"), min.cutoff = "q9")
FeaturePlot(vagus.combined, features = c("Phox2b"), min.cutoff = "q9")

### find all markers of cluster 3:

cluster3.markers <- FindMarkers(vagus.combined, ident.1 = 3, min.pct = 0.25)

### it lasts too long in my 8 Gb PC so I will switch to the PC with 
### higher RAM (96 Gb) .

cluster33.markers <- FindMarkers(vagus.combined, ident.1 = 33, min.pct = 0.25)


cluster30.markers <- FindMarkers(vagus.combined, ident.1 = 30, min.pct = 0.25)


cluster14.markers <- FindMarkers(vagus.combined, ident.1 = 14, min.pct = 0.25)


cluster2.markers <- FindMarkers(vagus.combined, ident.1 = 2, min.pct = 0.25)


cluster3.markers <- FindMarkers(vagus.combined, ident.1 = 3, min.pct = 0.25)

cluster7.markers <- FindMarkers(vagus.combined, ident.1 = 7, min.pct = 0.25)

cluster4.markers <- FindMarkers(vagus.combined, ident.1 = 4, min.pct = 0.25)


saveRDS(vagus.combined, file = "vagus_combine.rds")


write.csv(cluster14.markers, file = "cluster14.markers.csv")
write.csv(cluster2.markers, file = "cluster2.markers.csv")
write.csv(cluster3.markers, file = "cluster3.markers.csv")
write.csv(cluster7.markers, file = "cluster7.markers.csv")
write.csv(cluster4.markers, file = "cluster4.markers.csv")
write.csv(cluster30.markers, file = "cluster30.markers.csv")
write.csv(cluster33.markers, file = "cluster33.markers.csv") 

sessionInfo()
