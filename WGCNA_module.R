
## Performs module analysis on selected module
## Exports module eigengene expression
## Exports hub gene and eigengene lists
## Exports MM vs. Intramodular connectivity

## Requires x11

## ========================================================================
## Parameters to change
## ========================================================================
module = "turquoise"
datMEs = moduleEigengenes(datExpr_t, moduleColors)$eigengenes


## ========================================================================
## Heatmaps for Module Expression
## ========================================================================

quartz(width=10,height=5,pointsize=10)

which.module=module
ME=datMEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr_t[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2, cex.names=0.75,
        ylab="eigengene expression",xlab="array sample",names.arg=sample_times[order(sample_times_ordered_med)],axisnames=TRUE)
#axis(side=1, at=c(seq(.75,26.5,1.2)), labels=sample_names, tick=FALSE, padj=2,cex.axis=0.75)

quartz.save(file=paste("eigengene_expression_", softPower, "_merge", MEDissThres, "_module=", module, ".png", sep=""), type="png")

dev.off()

## ========================================================================
## Calculate Intramodular connectivity
## ========================================================================

ADJ1=abs(cor(datExpr_t,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)

## Add module color information to the Alldegrees table
Alldegrees1$moduleColors<-moduleColors

## Rank genes based on connectivity
color<-module
HubValues<-Alldegrees1[Alldegrees1$moduleColors==color,]
HubValues_ordered<-HubValues[order(HubValues$kWithin, decreasing=TRUE),]

HubValues_ordered_annotated = merge(HubValues_ordered[order(row.names(HubValues_ordered)),], blast_annotations, by.x=0, by.y=1, all.x=TRUE)

HubValues_ordered_annotated<-HubValues_ordered_annotated[order(HubValues_ordered_annotated$kWithin, decreasing=TRUE),]

write.table(HubValues_ordered_annotated, file=paste("hub_genes_beta", softPower, "_merge", MEDissThres, "_module=", module, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

## ========================================================================
## Locate genes closet to eigengene (Module Membership)
## ========================================================================

## Calculate module membership
datKME=signedKME(datExpr_t, datMEs, outputColumnName="MM.")

## Filter genes based on high module membership in a given module
mm_name = paste("MM.",module,sep="")
FilterGenes = abs(datKME[mm_name])>.8
table(FilterGenes) # How many were filtered out?
ModuleMembershipGeneList<-dimnames(data.frame(datExpr_t))[[2]][FilterGenes]
MembershipValue<-datKME[FilterGenes,mm_name]
MMTable<-cbind(ModuleMembershipGeneList,MembershipValue)
MMTableOrdered<-MMTable[order(MembershipValue,decreasing=TRUE),]
MMTableOrdered[,1]<-gsub("\\.", "|", MMTableOrdered[,1])

# Only get results that are within this module
MM_ordered = as.data.frame(MMTableOrdered)
MM_ordered = MM_ordered[MM_ordered$ModuleMembershipGeneList %in% all_info[mergedColors==module,1],]
MMTableOrdered = MM_ordered

MM_ordered_annotated = merge(MM_ordered[order(MM_ordered$ModuleMembershipGeneList),], blast_annotations, by.x=1, by.y=1, all.x=TRUE)
MM_ordered_annotated<-MM_ordered_annotated[order(MM_ordered_annotated$MembershipValue,decreasing=TRUE),]

write.table(MM_ordered_annotated, file=paste("genes_close_to_eigengenes_beta", softPower, "_merge", MEDissThres, "_module=", module, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

## ========================================================================
## Check relationship between MM and Intramdular Connectivity
## ========================================================================
x11(type="cairo")

which.color=module;
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

savePlot(filename=paste("intramodconn_beta", softPower, "_merge", MEDissThres, "_module=", module, ".png", sep=""), type="png")

dev.off()

## ========================================================================
## Plot module memberships and hub gene connectivities
## ========================================================================
# Is there a large drop-off in connectivity/MM?

quartz()

Hubs_df = cbind(Row.names=rownames(HubValues_ordered), HubValues_ordered)
Hubs_df = as.data.frame(cbind(Hubs_df$Row.names, Hubs_df$kWithin))
plot(Hubs_df[2], main="Hub Genes", xlab="Genes", ylab="Connectivity")

plot(MM_ordered[2], main="Genes closest to eigengenes", xlab="Genes", ylab="Module Membership")



## ========================================================================
## Get module membership values of ALL genes in a module with annotations
## ========================================================================

mod_mem = datKME[mm_name]
mod_mem = cbind(ids = rownames(mod_mem), mod_mem)
rownames(mod_mem) = NULL
mod_mem[,1]<-gsub("\\.", "|", mod_mem[,1])
mod_mem = mod_mem[mod_mem$ids %in% all_info[mergedColors==module,1],]

mod_mem = merge(mod_mem, blast_annotations, by.x=1, by.y=1, all.x=TRUE)
mod_mem = mod_mem[order(mod_mem[,2,drop=FALSE],decreasing=TRUE),]
rownames(mod_mem) = seq(1, nrow(mod_mem),1)

write.table(mod_mem, file=paste("mod_mem_b", softPower, "_merge", MEDissThres, "_module=", module, ".txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")


## Relationship between genes of interest and hub/genes close to eigengenes

# Rename row names from 1, 2, 3...
Hub_relation = HubValues_ordered_annotated 
rownames(Hub_relation) = seq(1,nrow(Hub_relation),1)

# The row number is the rank number
Hub_relation = Hub_relation[Hub_relation$Row.names %in% canonical_genes$V1,]
MM_relation = mod_mem[mod_mem$ids %in% canonical_genes$V1,]
