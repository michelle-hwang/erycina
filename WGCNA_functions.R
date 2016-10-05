##########################################################################################

get_adjacency_sums <- function(adj) {
  adjacency_sums = as.data.frame(rowSums(adj))
  adjacency_sums = cbind(Row.names=rownames(adjacency_sums), adjacency_sums)
  rownames(adjacency_sums) = NULL
  colnames(adjacency_sums) = c("Row.names", "sum_adjacency")
  return(adjacency_sums[order(adjacency_sums$sum_adjacency, decreasing=TRUE),])
  
}

##########################################################################################

get_one_module_variance <- function(MEDissThres) {
  merge = mergeCloseModules(datExpr_t, dynamicColors, cutHeight=MEDissThres, verbose=3)
  
  var_vals = as.data.frame(propVarExplained(datExpr_t, merge$colors, merge$newMEs))
  colnames(var_vals) = MEDissThres
  return(var_vals)
}

##########################################################################################

get_module_variances <- function(thresholds) {
  i = 0
  for(t in thresholds) {
    meow = get_one_module_variance(t)
    
    if(i==0) {
      meows = meow
    }
    else {
      meows = merge(meows, meow, by=0, all.x=TRUE, all.y=FALSE)
      meows2 = meows[,-1,drop=FALSE]
      rownames(meows2) = meows[,1]
      meows = meows2
    }
    
    i=i+1
  }
  return(meows2)
}

##########################################################################################

get_ranked_adjacencies <- function(gene) {
  adj = t(adjacency[gene,,drop=FALSE])
  adj = as.data.frame(adj[order(adj, decreasing=TRUE),,drop=FALSE])
  
  rownames(adj)<-gsub("\\.", "|", rownames(adj)) 
  adj_list = as.data.frame(adj)
  
  top_adj = merge(adj_list, all_info, by.x=0, by.y=1, all.x=TRUE, all.y=FALSE)
  top_adj = merge(top_adj, blast_annotations, by=1, all.x=TRUE)
  top_adj = top_adj[!top_adj$Row.names==gene,] # Remove itself
  top_adj[,2] = as.numeric(top_adj[,2])
  top_adj = top_adj[order(top_adj[,2], decreasing=TRUE),]
  
  return(top_adj)
}

##########################################################################################

get_ranked_adjacencies_overall <- function(genes, module) {
  # Build data frame
  genes_df = data.frame(matrix(ncol=0,nrow=5000))
  for( gene in genes ) {
    temp = t(adjacency[gene,,drop=FALSE])
    temp = as.data.frame(temp)
    genes_df = cbind(genes_df,temp)
  }
  rownames(genes_df)<-gsub("\\.", "|", rownames(genes_df))
  head(genes_df, n=15)
  
  # Get module probes
  probes = colnames(datExpr_t)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  
  # Use probe to remvoe unwanted adjacencies
  genes_df = genes_df[rownames(genes_df) %in% modProbes,]
  
  # Add col for row sums and row means
  genes_df = cbind(data.matrix(genes_df), rowSums(genes_df), rowMeans(genes_df))
  genes_df = as.data.frame(genes_df)
  colnames(genes_df) = c(genes, "sums", "means")
  
  # Add annotation
  genes_df = merge(genes_df, blast_annotations, by.x=0, by.y=1, all.x=TRUE)
  
  # Change column accordingly
  genes_df = genes_df[order(genes_df$sums, decreasing=TRUE),]
  head(genes_df, n=30)
  
  return(genes_df)
}

##########################################################################################

goi_rankings_in_hub_gene_lists <- function(module) {
  HubValues<-Alldegrees1[Alldegrees1$moduleColors==module,]
  HubValues_ordered<-HubValues[order(HubValues$kWithin, decreasing=TRUE),]
  
  df = merge(subset(all_results, mergedColors==module), HubValues_ordered, by.x=1, by.y=0, all.x=TRUE, all.y=FALSE)
  df = df[,c("colnames.datExpr_t.", "gene_set[annotations, 2]","kWithin")]
  colnames(df) = c("ids", "annotation", "kWithin")
  return(df)
}

##########################################################################################

module_trait_heatmap <- function() {
  MEs0<- moduleEigengenes(datExpr_t, moduleColors)$eigengenes
  MEs0_ordered<-orderMEs(MEs0)
  moduleTraitCor<-cor(MEs0_ordered, traitDat, use="p")
  moduleTraitPvalue<-corPvalueStudent(moduleTraitCor, nSamples)
  
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(traitDat),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}

##########################################################################################

plot_GSvsMM <- function(module, trait) {
  # Define variable weight containing the variable column of traitData
  trait_of_interest = as.data.frame(trait)
  
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr_t, MEs, use = "p"));
  
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(datExpr_t, trait_of_interest, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  names(geneTraitSignificance) = paste("GS.", names(trait_of_interest), sep="");
  names(GSPvalue) = paste("p.GS.", names(trait_of_interest), sep="");
  
  # Create plot of GS for variable vs MM for a given module and trait
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for trait",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}

##########################################################################################

plot_top_isoforms <- function(isoforms) {
  top = head(isoforms$V1, n=5)
  top = sub('\\Q|\\E', '\\.', top10,perl=TRUE)
  
  top_isoforms = NULL
  adjacency_sum_ordered = get_adjacency_sums(adjacency)
  for( x in top) {
    top_isoforms = rbind(top_isoforms,adjacency_sum_ordered[grep(x, adjacency_sum_ordered$Row.names),])
  }
  
  # Split up isoform ID to gene + isoform
  top_isoforms = cbind(data.frame(do.call('rbind', strsplit(as.character(top_isoforms$Row.names),'|',fixed=TRUE))), top_isoforms$sum_adjacency, top_isoforms$Row.names)
  
  # Simplify isoform ID to two chars
  top_isoforms$X2 = lapply(top_isoforms$X2, function(x) { substr(x,nchar(x)-1,nchar(x))})
  colnames(top_isoforms) = c("gene","isoform","sum_adjacency", "ids")
  
  quartz()
  
  (ggplot(top_isoforms, aes(x=ids, y=sum_adjacency, colour=factor(gene))) 
  + geom_point() 
  + ggtitle("Comparison of adjacency sums between isoforms") 
  + xlab("Isoforms") 
  + ylab("Adjacency Sums") 
  + scale_x_discrete(labels=top_isoforms$isoform)
  + scale_y_continuous(breaks=seq(0,1000,100))
  + theme(plot.title=element_text(size=15, vjust=1, face="bold"))
  + theme(axis.title.x=element_text(size=15, vjust=-0.5,hjust=0.5))
  + theme(axis.title.y=element_text(size=15, vjust=1.5,hjust=0.5))
  + theme(plot.margin=unit(c(1,1,1,1),"cm"))
  )
}
