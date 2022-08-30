library(dplyr)
library(Seurat)
library(patchwork)
dyn.load('/gpfs/data01/glasslab/home/zhl022/download/hdf5/lib/libhdf5_hl.so.200')
library(hdf5r)
library(harmony)

pbmc.data<-Read10X_h5("/data/scratch/zhl022/fetalMicroglia/all_sample_deep_Jun18_aggr/aggr_Jun18_2020/aggr/outs/raw_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = pbmc.data,  min.cells = 3, min.features = 200)
################## add sample name
#####################################################
pbmc@meta.data$sample_idx<-unlist(lapply(colnames(pbmc),function(x) unlist(strsplit(x,split="-"))[2] ))
con <- function(x){
  if (x=="1") {
    pre="A118"
  } else if (x=="2") {
    pre="42NC"
  } else if (x=="3") {
    pre="43NC"
  } else if (x=="4") {
    pre="44NC"
  } else if (x=="5") {
    pre="45NC"
  } else if (x=="6") {
    pre="46NC"
  } else if (x=="7") {
    pre="P89"
  } else if (x=="8") {
    pre="P90"
  } 
  return(pre)
}
pbmc@meta.data$sample<-unlist(lapply(pbmc@meta.data$sample_idx, con))
#####################################################
#####################################################
#####################################################
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca",group.by = "sample")

#pbmc <- FindNeighbors(pbmc, dims = 1:30)
#pbmc <- RunUMAP(pbmc, dims = 1:30)
#pbmc <- RunTSNE(pbmc, dims = 1:30)
#pbmc <- FindClusters(pbmc, resolution = 0.1)
                                         
pbmc <- RunHarmony(pbmc, "sample")
pbmc <- FindNeighbors(pbmc, reduction = "harmony",dims = 1:50)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, reduction = "harmony",dims = 1:50)
