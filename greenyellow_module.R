library(tidyverse)
library(WGCNA)

# first find genes in blue module
greenyellow <- gene_color_df[gene_color_df$moduleColors=="greenyellow",][,"genes"]
# subset combined_df by blue
greenyellow_df <- corrected_df[,greenyellow]
# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency = adjacency(greenyellow_df, power = softPower)

# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

colnames(TOM) <- greenyellow
rownames(TOM) <- greenyellow
colnames(dissTOM) <- greenyellow
rownames(dissTOM) <- greenyellow

# now choose the top connected genes for the genes of interest+ the top hub
greenyellow_genes_2_explore <- c('KDM4B','POLL')
# finding the nearest neighbors of a gene is to look up the top items in sorted columns of adjacency
# matrix

greenyellow_adjacency <- as.data.frame(adjacency[,greenyellow_genes_2_explore])
# stard8-first_gene
kdm4b <- greenyellow_adjacency[order(greenyellow_adjacency$KDM4B,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polb
poll <- greenyellow_adjacency[order(greenyellow_adjacency$POLL,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames

final_greenyellow <- c(kdm4b,poll)

# now to export module to cytoscape
modTOM = 1 - dissTOM
new_modTOM <- modTOM[final_greenyellow,final_greenyellow]
dimnames(new_modTOM) = list(final_greenyellow,final_greenyellow)
efname = "greenyellow_edges.txt"
nfname = "greenyellow_nodes.txt"

exportNetworkToCytoscape(new_modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = final_greenyellow,
)


