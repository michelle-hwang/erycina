
## Assumes "WGCNA_module.R" has been run

## ========================================================================
## Identifing important genes from DE between treatments
## ========================================================================
# At what rank are the canonical genes located?
# At what rank are the hub genes and genes closest to eigengenes located?
# Are they significantly differentially expressed?

# Note that this code will take ALL of your genes of interest, irrespective of module

wheres_waldo = function( DE_results ) {
  DE_results = cbind(ids=rownames(DE_results), DE_results)
  rownames(DE_results) = NULL
  DE_results$ids = sub('\\Q|\\E', '\\.', DE_results$ids, perl=TRUE)
  
  # Get hub genes
  Hub_ordered = cbind(ids=rownames(HubValues_ordered), HubValues_ordered)
  rownames(Hub_ordered) = NULL
  hub = cbind(head(Hub_ordered[1]), "hub")
  colnames(hub)[2] = "group"

  # Get genes closest to eigengenes
  M_ordered = MM_ordered
  eig = cbind(head(M_ordered[1]), "eigengene")
  colnames(eig) = c("ids", "group")
  
  # Get CAM/clock genes
  can = cbind(canonical_genes[1], "canonical")
  colnames(can) = c("ids", "group")
  
  # Put em together
  all = rbind(hub, eig, can)
  all$ids = sub('\\Q|\\E', '\\.', all$ids, perl=TRUE)

  # Now find the matches!
  matches = NULL
  for(x in all$ids) {
    matches = rbind(matches,DE_results[grep(x, DE_results$ids),])
  }
  return(merge(matches, all))
}


# Take matches and only return those that meet a p-value threshold
q <- function(de_results, pval) {
  results = wheres_waldo(de_results)
  results$PValue = sapply(results$PValue, as.numeric)
  results = results[results$PValue<pval,]
  
  # If no hits are found with p-value cut-off, return null
  if( nrow(results) == 0 ) {
    return(NULL)
  }
  
  # Add new col that specifies which pairwise comparison it is from
  x = cbind(results, deparse(substitute(de_results)))
  colnames(x)[7] = "pair"
  return(x)
}


# DE_results = output from edgeR
DE_2AM_vs_2PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.2AM_vs_2PM.edgeR.DE_results", header=TRUE)
DE_2AM_vs_6AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.2AM_vs_6AM.edgeR.DE_results", header=TRUE)
DE_2AM_vs_6PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.2AM_vs_6PM.edgeR.DE_results", header=TRUE)
DE_2PM_vs_6AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.2PM_vs_6AM.edgeR.DE_results", header=TRUE)
DE_2PM_vs_6PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.2PM_vs_6PM.edgeR.DE_results", header=TRUE)
DE_6AM_vs_6PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.6AM_vs_6PM.edgeR.DE_results", header=TRUE)
DE_10AM_vs_2AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10AM_vs_2AM.edgeR.DE_results", header=TRUE)
DE_10AM_vs_2PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10AM_vs_2PM.edgeR.DE_results", header=TRUE)
DE_10AM_vs_6AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10AM_vs_6AM.edgeR.DE_results", header=TRUE)
DE_10AM_vs_6PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10AM_vs_6PM.edgeR.DE_results", header=TRUE)
DE_10AM_vs_10PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10AM_vs_10PM.edgeR.DE_results", header=TRUE)
DE_10PM_vs_2AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10PM_vs_2AM.edgeR.DE_results", header=TRUE)
DE_10PM_vs_2PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10PM_vs_2PM.edgeR.DE_results", header=TRUE)
DE_10PM_vs_6AM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10PM_vs_6AM.edgeR.DE_results", header=TRUE)
DE_10PM_vs_6PM <- read.delim("../heatmaps/top5000/Erycina_top5000_trans.counts.matrix.10PM_vs_6PM.edgeR.DE_results", header=TRUE)


p = 0.05
significant_hits=NULL
significant_hits = rbind( q(DE_2AM_vs_2PM, p), 
                          q(DE_2AM_vs_6AM, p),
                          q(DE_2AM_vs_6PM, p),
                          q(DE_2PM_vs_6AM, p),
                          q(DE_2PM_vs_6PM, p),
                          q(DE_6AM_vs_6PM, p), 
                          q(DE_10AM_vs_2AM, p),
                          q(DE_10AM_vs_2PM, p),
                          q(DE_10AM_vs_6AM, p),
                          q(DE_10AM_vs_6PM, p),
                          q(DE_10AM_vs_10PM, p),
                          q(DE_10PM_vs_2AM, p),
                          q(DE_10PM_vs_2PM, p),
                          q(DE_10PM_vs_6AM, p),
                          q(DE_10PM_vs_6PM, p)
)

significant_hits[,1]<-gsub("\\.", "|", significant_hits[,1])
significant_hits = merge(significant_hits, blast_annotations, by=1, all.x=TRUE)

write.table(significant_hits, file="significant_hits-antiquewhite4.txt", sep="\t", quote=FALSE, row.names=TRUE)




