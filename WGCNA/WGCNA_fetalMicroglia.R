library(WGCNA)
library(flashClust)
library(dplyr)
library('AnnotationDbi')
library("org.Hs.eg.db")
library(ggplot2)
library(igraph)




datExpr = read.csv("20200821_Master_TPM.txt",sep='\t',row.names = 1)



datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
dim(datExpr) # 48 samples and 1000 genes (you will have many more genes in reality)


# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK


#If the last statement returns TRUE, all genes have passed the cuts. 
#If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
  {if (sum(!gsg$goodGenes)>0)
       printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
       if (sum(!gsg$goodSamples)>0)
           printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
       datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
       }


#Create an object called "datTraits" that contains your trait data
datTraits = read.csv("trait_WGCNA_1.csv" )
datTraits<- datTraits %>% mutate(genders = ifelse(gender %in% c("Male"),1,0))
rownames(datTraits) = datTraits$sample
datTraits$sample = NULL
datTraits$gender = NULL
head(datTraits)
datExpr <- datExpr[rownames(datExpr) %in% rownames(datTraits),]
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)

###########################################################################Select Soft Power

# Choose a soft threshold power- USE A SUPERCOMPUTER IRL ------------------------------------

powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

#from this plot, we would choose a power of 18 becuase it's the lowest power for which the scale free topology index reaches 0.90

###########################################################################Construct a gene co-expression matrix and generate modules
#build a adjacency "correlation" matrix
enableWGCNAThreads()
softPower = 10
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
head(adjacency)

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# Generate Modules --------------------------------------------------------


# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEs[, colSums(is.nan(MEs)) != nrow(MEs)]
MEs_rm_nan<-MEs %>% select_if(~ !any(is.nan(.)))
MEs<-MEs_rm_nan

MEDiss= 1-cor(MEs, use="pairwise.complete.obs")
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")


#plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="clusterwithoutmodulecolors.pdf")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.01
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

# Correlate traits --------------------------------------------------------


#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))


#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

names(datExpr)[moduleColors=="brown"]

modules = c("darkorange");

probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

modTOM = TOM[inModule, inModule];

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE, threshold = 0,
                               nodeNames = modProbes, nodeAttr = moduleColors[inModule])

rownames(modTOM) <- modProbes
colnames(modTOM) <- modProbes
########################################### compute hub score
gra<-igraph::graph_from_adjacency_matrix(modTOM, mode = c( "undirected"),weighted = T)
hub_genes<-hub_score(gra, scale = TRUE)

