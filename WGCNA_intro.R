
## =============================================================================
# WGCNA Analysis
## =============================================================================

## Full set of scripts used to run a WGCNA analysis given a counts matrix (for Pearson correlation) or a similarity/adjacency matrix calculated from rsgcc package (for gini correlation) on a time series experiment.

## WGCNA_intro.R
## WGCNA_analysis.R
## WGCNA_module.R
## WGCNA_deresults.R
## WGCNA_cytoscape.R
## WGCNA_misc.R

## Adapted by Michelle Hwang and Amanda Kelly Lane
## Original code from "Tutorials for the WGCNA package" by Peter Langfelder and Steven Horvath

# See WGCNA documentation to install WGCNA package
library(WGCNA)
library(edgeR)
library(rsgcc)

options(stringsAsFactors=FALSE)
setwd("~/Google Drive/UGA/Research/Erycina\ 2/")

## Used for plot labels
sample_names <- c('1','3','4','5','6','7','8','9','11','13','15','18','25','26','27','28','29','30','31','32','33','34','35','36')
sample_times <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '2AM', '4PM', '10PM','8PM', '6PM', '6AM', '12PM', '4AM', '4PM', '12PM', '2AM', '8AM', '2PM', '10AM', '8AM', '2PM', '6PM', '8PM', '12AM'))
sample_times_v2 <- c(22, 0, 4, 6, 10, 2, 16, 22, 20, 18, 6, 12, 4, 16, 12, 2, 8, 14, 10, 8, 14, 18, 20, 0)
cbind(sample_names[order(sample_times_v2)], sample_times[order(sample_times_v2)])
#sample_names <- c('1','3','4','5','6','7','8','9','11','13','15','18','25','26','27','28','29','30','31','32','33','34','35','36', 'E3','E4','E5','E6','E7','E8','E9','E10','E12','E14','E15','E16','E17','E18','E19','E20','E21','E43','E44','E45','E59','E60')
#sample_times <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '2AM', '4PM', '10PM','8PM', '6PM', '6AM', '12PM', '4AM', '4PM', '12PM', '2AM', '8AM', '2PM', '10AM', '8AM', '2PM', '6PM', '8PM', '12AM'), c('2AM', '2AM','6PM','10PM','2PM','2AM','2PM','2PM','2PM','6PM','6PM','6PM','6PM','10PM','10PM','10PM','10PM','10AM','10AM','10AM','6AM','6AM'))
#sample_times_v2 <- c(c(22, 0, 4, 6, 10, 2, 16, 22, 20, 18, 6, 12, 4, 16, 12, 2, 8, 14, 10, 8, 14, 18, 20, 0), c(2, 2, 18, 22, 14, 2, 14, 14, 14, 18, 18, 18, 18, 22, 22, 22, 22, 10, 10, 10, 6, 6))
#cbind(sample_names[order(sample_times_v2)], sample_times[order(sample_times_v2)])

sample_names <- c('1','3','4','5','7','8','9','13','15','18','25','26','29','31','32','33','34','35','36', 'E3','E4','E6','E7','E8','E9','E10','E12','E14','E15','E16','E17','E18','E19','E20','E21','E43','E44','E45','E59','E60')
sample_times <- c(c('10PM', '12AM', '4AM', '6AM', '10AM', '4PM', '10PM','8PM', '6PM', '12PM', '4AM', '4PM',  '8AM', '10AM', '8AM', '2PM', '6PM', '8PM', '12AM'), c('2AM', '2AM','10PM','2PM','2AM','2PM','2PM','2PM','6PM','6PM','6PM','6PM','10PM','10PM','10PM','10PM','10AM','10AM','10AM','6AM','6AM'))
sample_times_v2 <- c(c(22, 0, 4, 6, 10, 16, 22, 20, 18, 12, 4, 16, 20, 22, 8, 14, 10, 8, 14, 18, 20, 0), c(2, 2, 22, 14, 2, 14, 14, 14, 18, 18, 18, 18, 22, 22, 22, 22, 10, 10, 10, 6, 6))
cbind(sample_names[order(sample_times_v2)], sample_times[order(sample_times_v2)])

## ========================================================================
## Get all the data you need
## ========================================================================

## Import counts matrix: rows=samples, cols=genes
datExpr = read.delim("data/feb25/d_ndm_all_5000_newcorr_with_i_top5000.counts.matrix", row.names=1)
datExpr = read.delim("data/feb25/d_ndm_10000_newcorr.counts.matrix", row.names=1)
colnames(datExpr) = sample_names
datExpr = datExpr[,order(sample_times_ordered_med)]
datExpr_t <- t(datExpr) # Transpose if needed
nGenes = ncol(datExpr_t)
nSamples = nrow(datExpr_t)
datExprTMM <- read.delim("data/feb25/d_ndm_all_5000_newcorr_with_i_top5000.txt", row.names=1)
datExprTMM <- read.delim("data/feb25/d_ndm_all_10000_newcorr.txt", row.names=1)

datExprG = adjacencymatrix(datExpr, method="GCC")


## Import adjacency/similarity matrix (gini)
#corMat <- read.delim("data/feb9/diffExpr.P0.001.iso4.matrix", row.names=1, stringsAsFactors=FALSE)

## Import trait data: rows=samples, cols=traits
#traitDat <- read.delim("data/trait_data.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

## List of genes of interest: col1=gene_id, col2=gene_annotation
# Tab-delimited format: ID + Annotation
canonical_genes <- read.delim('data/canonical_genes.txt', header=FALSE)
cam_genes <- read.delim('data/cam_genes.txt', header=FALSE)

## BLAST Annotations
# Tab-delimited format: ID + Annotation
blast_annotations <- read.delim('data/batch_blastx_annotations_short.txt', header=FALSE)
blast_annotations = unique(blast_annotations)
blast_annotations = blast_annotations[order(blast_annotations$V1),]
colnames(blast_annotations) = c("ID", "annotation")


## ========================================================================
## Dendrogram of Samples
## ========================================================================
quartz()
tree = hclust(dist(t(datExprTMM)), method="average")
tree$labels = sample_times
col_maps = as.data.frame(cbind(c(rep("darkblue",24),sample_times)))
col_maps = col_maps[tree$order,]
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    i <<- i+1
    attr(n, "nodePar") <- list(a$nodePar, lab.col=col_maps[i,1], lab.cex=0.8)
    print(i)
  }
  n
}
i=0
tree2 <- dendrapply(as.dendrogram(tree), colLab)
plot(tree, main="cluster dendrogram of samples", horiz=TRUE)
plot(tree, main="cluster dendrogram of samples")
#quartz.save(file="cluster_dendrogram_of_samples.png")
## ========================================================================
## Choosing soft thresholding power
## ========================================================================
choose_beta <- function(cor_mat) {
  powers = c(seq(from=1, to=20, by=1))
  #sft = pickSoftThreshold(datExpr_t, powerVector=powers, verbose=5, RsquaredCut=0.7)
  sft = pickSoftThreshold.fromSimilarity(datExprG, powerVector=powers, verbose=5, RsquaredCut=0.7,  networkType = "unsigned")
  
  x11(type="cairo")
  par(mfrow=c(1,2))
  cex1=0.9
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",
       ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}


choose_beta(corMat)


## MDS plot 
plotMDS(datExprTMM, labels=paste(sample_times,sample_names,sep="-"))
plotMDS(datExprTMM, col=c(rep("black",24),rep("red",22)))
plotMDS(datExprTMM, col=c(rep("black",24),rep("red",22)), pch=20, labels)

plotMDS(datExprTMM, labels=paste(sample_times, sep="-"))
plotMDS(datExprTMM, col=col_maps[,1], pch=20, labels=sample_times_v2)
quartz()




