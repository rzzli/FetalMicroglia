library(velocyto.R)
dyn.load('/gpfs/data01/glasslab/home/zhl022/download/hdf5/lib/libhdf5_hl.so.200')
library(pagoda2)
library(Seurat)

ldat <- read.loom.matrices("merge_Jun21.loom") #loom file is available upon request


DimPlot(pbmc,reduction = 'tsne')

###################################
## change names
#

cellnames<-colnames(ldat$spliced) # get cell names from loom
cell_pre <-unlist(lapply(cellnames, function(x) (unlist(strsplit(x,split=':'))[1])))
cell_pre_unique <- unique(cell_pre)
print(cell_pre_unique) #loom pre
####
#
cluster_tsne <- as.data.frame(fetal_mg_sc_seurat@reductions$tsne@cell.embeddings)
cluster_tsne$cluster <-  as.character(fetal_mg_sc_seurat@meta.data$seurat_clusters)

#con convert names from Seurat to Velocyto
con <- function(rn){
  bar<-unlist(strsplit(rn,split = "-"))[1]
  pre_sample<-unlist(strsplit(rn,split = "-"))[2]
  if (pre_sample=="1") {
    pre="A118:"
  } else if (pre_sample=="2") {
    pre="JB_134_14w1d_usp_42NC:"
  } else if (pre_sample=="3") {
    pre="JB_136_12w1d_ctx_44NC:"
  } else if (pre_sample=="4") {
    pre="JB_207:"
  } else if (pre_sample=="5") {
    pre="JB_135_12w1d_bs_43NC:"
  } else if (pre_sample=="6") {
    pre="JB_138_10w1d_usp_45NC:"
  } else if (pre_sample=="7") {
    pre="JB_139_12w6d_usp_46NC:"
  } else if (pre_sample=="8") {
    pre="JB_235:"
  } 
  new <-paste(pre,bar,'x',sep="")
  return(new)
}

seurat_to_loom_cell_name <- unlist(lapply(rownames(cluster_tsne),con))
rownames(cluster_tsne)<-seurat_to_loom_cell_name
###########################################################
###########################################################

# this dataset has already been pre-filtered, but this is where one woudl do some filtering
emat <- ldat$spliced

emat <- emat[,colSums(emat)>=3e3]
celln<-dim(emat)[2]
#cellnr<-sample(1:celln,50000)
#emat <- emat[,cellnr]

rownames(emat) <- make.unique(rownames(emat))
########################################################################################
##
#reassign cluster and reassign tsne
cells <- colnames(emat)
tsne <- cluster_tsne ## reassign name to tsne
tsne_sub<-subset(tsne, rownames(tsne) %in% cells)
##
########################################################################################



##re_subset  emat
emat<-emat[,rownames(tsne_sub)]


r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

r$adjustVariance(plot=T,do.par=T,gam.k=10)

r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)

r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')


r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
########################################################################################
##
test_cluster<-r$clusters$PCA[[1]]
test_cluster_x<-as.factor(tsne_sub[,3])
names(test_cluster_x)<-rownames(tsne_sub) # cluster assignment
# now test_cluster_x is the cluster annotation

test_tsne<-as.matrix(tsne_sub[,1:2])  #tsne 
# now test_tsne is the official coordinates
#r$embeddings$PCA$tSNE
#r@.xData$embeddings$PCA$tSNE
#r$embeddings$PCA$tSNE <- NULL
#r$clusters$PCA$multilevel
r$clusters$PCA$multilevel <-test_cluster_x
r$embeddings$PCA$tSNE<-test_tsne
##
########################################################################################

# back to plotting

par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"APOE"],main='APOE')




#####Velocity estimation
emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE

cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
par(mfrow=c(1,1))
pdf("velocyto_Oct20.pdf")
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()
