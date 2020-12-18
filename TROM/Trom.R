library(lattice)
library(gplots)
library(RColorBrewer)
library(gtools)
library(openxlsx)
library(topGO)
library(GO.db)
library(TROM)
library(gdata)
library(pheatmap)

mg <- read.table("FP_TPM_20201119.txt",header = T,row.names = 1)
rownames(ctx)  <- make.unique(as.character(ctx$X))
genes<-rownames(mg)
mg <-mg[,1:88]

mg <-log2(mg+1)

dm_trom<-ws.trom(sp_gene_expr = mg,provide = F,z_thre=0.75 , save_overlap_genes = F)
mat <- ifelse(dm_trom>6,6,dm_trom)

pdf("similarity_map_Dec7.pdf")
pheatmap(mat,cluster_cols = F,cluster_rows = F,fontsize =5)
dev.off()
