library(tidyverse)
library(WGCNA)
# first find genes in blue module
# first find genes in blue module
# first find genes in blue module
turquoise <- gene_color_df[gene_color_df$moduleColors=="turquoise",][,"genes"]
# subset combined_df by blue
turquoise_df <- corrected_df[,turquoise]
# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency = adjacency(turquoise_df, power = softPower)

# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

colnames(TOM) <- turquoise
rownames(TOM) <- turquoise
colnames(dissTOM) <- turquoise
rownames(dissTOM) <- turquoise

# now choose the top connected genes for the genes of interest+ the top hub
turquoise_genes_2_explore <- c('LRRC1','REV1')
# finding the nearest neighbors of a gene is to look up the top items in sorted columns of adjacency
# matrix

turquoise_adjacency <- as.data.frame(adjacency[,turquoise_genes_2_explore])
# lrrc1-first_gene
lrrc1 <- turquoise_adjacency[order(turquoise_adjacency$LRRC1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# rev1
rev1 <- turquoise_adjacency[order(turquoise_adjacency$REV1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames

final_turquoise <- c(lrrc1,rev1)

# now to export module to cytoscape
modTOM = 1 - dissTOM
new_modTOM <- modTOM[final_turquoise,final_turquoise]
dimnames(new_modTOM) = list(final_turquoise,final_turquoise)
efname = "turquoise_edges.txt"
nfname = "turquoise_nodes.txt"

exportNetworkToCytoscape(new_modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = final_turquoise,
)


