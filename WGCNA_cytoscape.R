
## Assumes "WGCNA_quick.R" has already been run.
## Runs instantaneously

## ========================================================================
## Parameters to change
## ========================================================================

## Select modules
## Current version only works with 1 module
modules = c("brown");
adjThresh = 0.5

## ========================================================================
## Main Code
## ========================================================================

## Get hub genes
ADJ1=abs(cor(datExpr_t,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
Alldegrees1$moduleColors<-moduleColors
HubValues<-Alldegrees1[Alldegrees1$moduleColors==modules,]
HubValues_ordered<-HubValues[order(HubValues$kWithin, decreasing=TRUE),]

## Get genes with highest module membership
datMEs = moduleEigengenes(datExpr_t, moduleColors)$eigengenes
datKME=signedKME(datExpr_t, datMEs, outputColumnName="MM.")
mm_name = paste("MM.",modules,sep="")
FilterGenes= abs(datKME[mm_name])>.8
ModuleMembershipGeneList<-dimnames(data.frame(datExpr_t))[[2]][FilterGenes]
MembershipValue<-datKME[FilterGenes,mm_name]
MMTable<-cbind(ModuleMembershipGeneList,MembershipValue)
MMTableOrdered<-MMTable[order(MembershipValue,decreasing=TRUE),]
MMTableOrdered[,1]<-gsub("\\.", "|", MMTableOrdered[,1])
MM_ordered = as.data.frame(MMTableOrdered)
MM_ordered = MM_ordered[MM_ordered$ModuleMembershipGeneList %in% all_info[mergedColors==module,1],]
MMTableOrdered = MM_ordered

## Get genes with highest connectivity to genes of interest
gois = subset(all_results, mergedColors==modules)[1]
gois = gois$colnames.datExpr_t.

unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

hub_connect = NULL
for(i in 1:length(gois)) {
  new_adjacency = unAsIs(adjacency)
  colnames(new_adjacency) = NULL
  row = new_adjacency[gois[i],]
  df = t(rbind(colnames(datExpr_t), row))
  df_ordered = df[order(row, decreasing=TRUE),]
  hub_connect = rbind(hub_connect, cbind( gois[i], head(n=10, df_ordered[,1]), "high"))
}
colnames(hub_connect) = c("fromNode","toNode","connection")

## Select module probes
probes = colnames(datExpr_t)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

## Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

## Annotate if CAM/clock, hub, or high module membership gene
annotProbes = c()
for( x in 1:length(modProbes)) {
  if( is.element(modProbes[x], cam_genes[,1])) {
    annotProbes = c(annotProbes, "1")
    if( is.element(modProbes[x], head(n=10, row.names(HubValues_ordered)))) {
      annotProbes = c(annotProbes, "4")
    }
    else if( is.element(modProbes[x], head(n=10, MMTableOrdered[,1]))) {
      annotProbes = c(annotProbes, "5")
    }
  }
  else if( is.element(modProbes[x], head(n=10, row.names(HubValues_ordered)))) {
    if( is.element(modProbes[x], head(n=10, MMTableOrdered[,1]))) {
      annotProbes = c(annotProbes, "6")
    }
    else {
      annotProbes = c(annotProbes, "2")
    }
  }
  else if( is.element(modProbes[x], head(n=10, MMTableOrdered[,1]))) {
    annotProbes = c(annotProbes, "3")
  }
  else {
    annotProbes = c(annotProbes, "NA")
  }
}

## Get the top x genes (ranked by adjacency)
#nTop = 1000;
#IMConn = softConnectivity(datExpr_t[, modProbes]);
#top = (rank(-IMConn) <= nTop)
#modTOM = modTOM[top,top]

modTOM[modTOM>1] <- 1

## Output Cytoscape import files
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Cyto-edges_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("Cyto-nodes_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = adjThresh,
                               nodeNames = modProbes,
                               altNodeNames = annotProbes,
                               nodeAttr = moduleColors[inModule]);

## Add edge annotations for edges with highest connectivity to CAM/clock genes
edgefile = read.delim(paste("Cyto-edges_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse="-"), ".txt", sep=""))

edgefile = merge(edgefile, hub_connect, by=c("fromNode","toNode"), all.x=TRUE)

write.table(edgefile, file=paste("Cyto-edges_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse='-'), ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

## Add BLASTX annotations to node annotations
nodefile = read.delim(paste("Cyto-nodes_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse="-"), ".txt", sep=""))

nodefile = merge(nodefile, blast_annotations, by=1, all.x=TRUE)

write.table(nodefile, file=paste("Cyto-nodes_b", paste(softPower, "_merge", MEDissThres, "_adj", adjThresh, "_", sep=""), paste(modules, collapse="-"), ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

