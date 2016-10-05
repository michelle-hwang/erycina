
## Assumes "WGCNA_analysis.R" has already been run.

library(grid)


## ========================================================================
## Module-Trait Analysis
## ========================================================================


## Correlate module-trait relationships
# Recalculate MEs with color labels
module_trait_heatmap()
dev.copy2pdf(file="module-acidity-relationships.pdf")


## Quantify associations of individual genes with the trait of interest
plot_GSvsMM("green", traitDat$acidity)
plot_GSvsMM("cyan", traitDat$acidity)
dev.copy2pdf(file="GSvsMM.pdf")


## Get the GSP for CAM/clock genes
geneTraitSignificance = as.data.frame(cor(datExpr_t, as.data.frame(traitDat$acidity), use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(GSPvalue) = "acidity"
goi_acidity = merge(all_results, GSPvalue, by.x=1, by.y=0, all.x=FALSE, all.y=FALSE)
goi_acidity = goi_acidity[order(goi_acidity$acidity),]

write.table(goi_acidity, file="gsp_goi.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Get table of GSP
write.table(GSPvalue[order(rownames(GSPvalue)),,drop=FALSE], file="gspvalues.txt", sep="\t", quote=FALSE, row.names=TRUE)


## ========================================================================
## Adjacency Analysis
## ========================================================================


## Get list of all genes ranked by adjacency sum 
adjacency_sum_ordered = get_adjacency_sums(adjacency)
write.table(adjacency_sum_ordered[order(adjacency_sum_ordered$Row.names),,drop=FALSE], file="adjacency_sum_ordered.txt", sep="\t", quote=FALSE, row.names=FALSE)

## Plot distribution of adjacency sums
quartz(width=9, height=4)

hist(adjacency_sum_ordered[2], breaks=60, main="Distribution of Adjacency Sums", xlab="Sum of adjacencies across isoforms", axes=FALSE, plot=TRUE)
axis(1, at=seq(0,1500,50))
axis(2, at=seq(0,2000,100))
quartz.save(file="adjacency_beta8.png", type="png")


## Get ranked list of adjacencies given gene, with annotations and module info
# Irrespective of module
gene="TR9881|c0_g1_i1"
list = get_ranked_adjacencies(gene)
head(n=20, list)
write.table(list, "adj_rank_TR9881|c0_g1_i1.txt", sep="\t", quote=FALSE, row.names=FALSE)


## Given 2+ CAM genes in a module, rank adjacency sums/means
genes = c("TR51771|c0_g1_i1","TR101322|c8_g1_i4","TR43986|c1_g1_i1","TR105212|c0_g1_i1","TR9881|c0_g1_i1")
module = "brown"

genes_df = get_ranked_adjacencies_overall(genes, module)
write.table(genes_df, file="top_adj_genes-brown.txt", quote=FALSE,row.names=FALSE, sep="\t")


## How are CAM genes ranked in hub gene lists? (connectivity lists)
module = "darkorange2"
goi_rankings_in_hub_gene_lists(module)


## Hetero rank vs. Adj sum 
id_df = read.delim("../data/top5000_DEresults_ids.txt", header=FALSE)
temp = as.data.frame(cbind(id_df, seq(1,5000,1)))
temp = merge(temp, adjacency_sum_ordered, by=1, all.x=TRUE)
temp$V2 = as.numeric(temp$V2)
rank_and_adj = temp[order(temp[,2]),]
plot(rank_and_adj[,2], rank_and_adj$sum_adjacency, main="DE Rank vs Adjacency", xlab="Heterogeneous expression rank", ylab="Adjacency sum")

## ========================================================================
## Isoform Analysis
## ========================================================================


## Format is tab-delimited: ID + Isoform count
isoforms_full = read.delim("../data/isoforms_full_set.txt", header=FALSE)
isoforms_5000 = read.delim("../data/isoforms_top5000.txt", header=FALSE)
isoforms_10000 = read.delim("../data/isoforms_top10000.txt", header=FALSE)

# Check out the data
table(isoforms_full$V2)
table(isoforms_5000$V2)
table(isoforms_10000$V2)


## Isoform distribution
isoform_data = cbind(as.vector(table(isoforms_full$V2)), as.vector(table(isoforms_10000$V2)), as.vector(table(isoforms_5000$V2))) # Combine the tables
isoform_data[6,3] = 0 # Add a 6,0 point to the last table (missing)
colnames(isoform_data) = c("full","10,000","5,000")
rownames(isoform_data) = seq(1,6,1)

barplot(isoform_data, beside=TRUE, axisnames=TRUE, main="Frequency of Isoform Counts", legend.text=TRUE, col=topo.colors(6))
quartz.save(file="plot_isoform_counts.png", type="png")


## Genes with high isoform counts with adjacency sums
plot_top_isoforms(isoforms_5000)
quartz.save(file="top_isoforms_adjacency_comparison.png", type="png")


## A gene's isoforms plotted against adjacency to all other genes
gene = "TR76619"
indices = grep(gene, rownames(adjacency))
cl=rainbow(length(indices)-1)
plot(adjacency[indices[1],], col=rgb(0,0,0,alpha=0.5), pch=19, main="Isoform Adjacencies of Gene TR76619", xlab="Genes", ylab="Adjacencies")
col_index = 1
for( index in indices[-1] ) {
  color = cl[col_index]
  color_trans = adjustcolor(color, alpha.f = 0.3)
  points(adjacency[index,], col=color_trans, pch=19)
  col_index = col_index + 1
}

quartz.save(file="isoform_adj_TR76619.png", type="png")


## ========================================================================
## Proportion of variance explained by eigengenes
## ========================================================================
thresholds = seq(0,0.45,0.025)
mod_vars = get_module_variances(thresholds)
mod_vars
mod_vars[is.na(mod_vars)] <- 0

quartz()
plot(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEturquoise",], col="turquoise", pch=19, main="Variance explained by eigenegenes at varying height cut-offs", xlab="Height Cut-off", ylab="Variance explained")
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEdarkorange2",], col="darkorange2", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEivory",], col="black", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEgreen",], col="green", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEantiquewhite4",], col="antiquewhite4", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEbrown",], col="brown", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEblue",], col="blue", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEdarkseagreen4",], col="darkseagreen4", pch=19)
points(x=as.numeric(colnames(mod_vars)), y=mod_vars["PVEsalmon4",], col="salmon4", pch=19)

colors = c("turquoise", "darkorange2", "black", "green", "antiquewhite4", "brown", "blue", "darkseagreen4")
legend( 0, 0.19, colors, col=colors, pch=19, pt.cex=1, cex=.85)

quartz.save("variance_explained.png")

## ========================================================================
## Other
## ========================================================================

## Custom output table
fpkm <- read.delim("../data/Erycina_trinity_trans_filtered.TMM.fpkm.matrix", row.names=1)

# Add an index - is it a canonical gene?
canonical_index = canonical_genes
canonical_index[,2] = rep(1, length(canonical_index[,2]))
FINAL_TABLE = merge(all_info, canonical_index, by=1, all.x=TRUE, all.y=FALSE)

# Add annotations
FINAL_TABLE = merge(FINAL_TABLE, blast_annotations, by=1, all.x=TRUE, all.y=FALSE)
colnames(FINAL_TABLE) = c("ids", "module_color", "module_num", "cam_or_clock", "annotation")

# Add fpkm counts
FINAL_TABLE = merge(FINAL_TABLE, datExpr, by.x=1, by.y=0, all.x=TRUE, all.y=FALSE)

# Add acidity
FINAL_TABLE = merge(FINAL_TABLE, GSPvalue[order(rownames(GSPvalue)),,drop=FALSE], by.x=1, by.y=0, all.x=TRUE)

write.table(FINAL_TABLE, "erycina_wgcna.txt", sep="\t", quote=FALSE, row.names=FALSE)






