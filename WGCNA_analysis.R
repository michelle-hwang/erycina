
## Quickly run WGCNA on a similarity or counts matrix. 
## Exports dendrogram + module before/after merge plot
## Exports dendrogram before/after merge cut-off
## Exports data table of all genes and their module assignments + annotations
## Exports data table of canonical genes and their module assignments + annotations
## Exports TOM plot 

## Runs 2-5 minutes dep. on size
## Requires x11

## ========================================================================
## Parameters to change
## ========================================================================
softPower = 9 #Beta thresholding power
minModuleSize = 20 #Min num genes to make up a module
MEDissThres = 0.02 #Post dynamic cut merge height threshold

## ========================================================================
## Co-expression similarity and adjacency
## ========================================================================
#adjacency = adjacency(datExpr_t, power=softPower)
adjacency = adjacency.fromSimilarity(datExprG, power=softPower, type="unsigned")

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method="average")

## ========================================================================
## Module identification using dynamic tree cut
## ========================================================================
dynamicMods = cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
nrow(table(dynamicColors))

## ========================================================================
## Merge modules w/ similar expression profilest
## ========================================================================
MEList = moduleEigengenes(datExpr_t, colors=dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs) # Calc dissimilarity of module eigengenes
METree = hclust(as.dist(MEDiss), method="average") # Cluster module eigengenes

x11(type="cairo")
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")
abline(h=MEDissThres, col="red")
dev.copy2pdf(file=paste("dendro_beta", softPower, "_merge_", MEDissThres, ".pdf", sep=""))

## Call automatic merging function
merge = mergeCloseModules(datExpr_t, dynamicColors, cutHeight=MEDissThres, verbose=3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
print("Num modules after merge:", quote=FALSE)
ncol(mergedMEs)

## Should any more modules be merged based on module membership?
moduleMergeUsingKME(datExpr_t, mergedColors)

## Save variables for later
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

## Make new dendrogram
MEDiss = 1-cor(MEs) # Calc dissimilarity of module eigengenes
METree = hclust(as.dist(MEDiss), method="average") # Cluster module eigengenes
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")
dev.copy2pdf(file=paste("dendro2_beta", softPower, "_merge_", MEDissThres, ".pdf", sep=""))

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)

savePlot(filename=paste("cluster_beta", softPower, "_merge", MEDissThres, ".png", sep=""), type="png")


## ========================================================================
## Exporting genes and networks
## ========================================================================

## Get dataframe with module info for entire data set
all_info <- data.frame(colnames(datExpr_t), mergedColors, moduleLabels)

## Function to subset data given gene list
get_genes <- function( gene_set, all_genes ){
  x=1
  gene_list <- NULL
  annotations <- NULL
  for( k in 1:nrow(all_info[1])) {
    for( i in 1:nrow(gene_set[1])) {
      if( gene_set[i,1] == all_genes[k,1]) {
        gene_list[x] <- k
        annotations[x] <- i
        x <- x+1
      }
    }
  }
  return( cbind(all_genes[gene_list,],gene_set[annotations,2]) )
}

## Rate limiting step
all_results <- get_genes(canonical_genes, all_info)
all_results <- unique(all_results)


## Output table for canonical genes
write.table(all_results, file=paste("modules_canonical_beta", softPower, "_merge", MEDissThres, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

## Output table for all genes
write.table(all_info, file=paste("modules_beta", softPower, "_merge", MEDissThres, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

## Output table for module counts
write.table(table(moduleColors), file=paste("module_counts_beta", softPower, "_merge", MEDissThres, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

## Visualizing the gene network
plotTOM = dissTOM
diag(plotTOM) = NA
geneTree = hclust(as.dist(dissTOM), method="average")

#TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes", col=c('red', 'orange','yellow','green','blue'), breaks=c(0,0.50,0.85,0.98,0.99,1.0))
#TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes", col=c('red', 'orange','yellow','green','blue'), breaks=c(0,0.60,0.70,0.80,0.95,1.0))

#savePlot(filename=paste("heatmap_beta", softPower, "_merge", MEDissThres, ".png", sep=""), type="png")

#dev.off()

## ========================================================================
  
## Plot pairwise scatter plot of samples along eigengene
## ========================================================================

#datMEs = moduleEigengenes(datExpr_t, moduleColors)$eigengenes

# Use this plot to check whether any of the modules are highly correlated with another
#quartz(4,10)
#plotMEpairs(datMEs, y=NULL)




