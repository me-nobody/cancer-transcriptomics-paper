library(tidyverse)
library(WGCNA)

# first find genes in blue module
red <- gene_color_df[gene_color_df$moduleColors=="red",][,"genes"]
# subset combined_df by blue
red_df <- corrected_df[,red]
# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency = adjacency(red_df, power = softPower)

# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

colnames(TOM) <- red
rownames(TOM) <- red
colnames(dissTOM) <- red
rownames(dissTOM) <- red

# now choose the top connected genes for the genes of interest+ the top hub
red_genes_2_explore <- c('VAX1','POLM')
# finding the nearest neighbors of a gene is to look up the top items in sorted columns of adjacency
# matrix

red_adjacency <- as.data.frame(adjacency[,red_genes_2_explore])
# stard8-first_gene
vax1 <- red_adjacency[order(red_adjacency$VAX1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polb
polm <- red_adjacency[order(red_adjacency$POLM,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames

final_red <- c(vax1,polm)

# now to export module to cytoscape
modTOM = 1 - dissTOM
new_modTOM <- modTOM[final_red,final_red]
dimnames(new_modTOM) = list(final_red,final_red)
efname = "red_edges.txt"
nfname = "red_nodes.txt"

exportNetworkToCytoscape(new_modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = final_red,
)


