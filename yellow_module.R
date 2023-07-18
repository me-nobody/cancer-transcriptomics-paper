library(tidyverse)
library(WGCNA)

# first find genes in blue module
yellow <- gene_color_df[gene_color_df$moduleColors=="yellow",][,"genes"]
# subset combined_df by blue
yellow_df <- corrected_df[,yellow]
# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency = as.data.frame(adjacency(yellow_df, power = softPower))

# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

colnames(TOM) <- yellow
rownames(TOM) <- yellow
colnames(dissTOM) <- yellow
rownames(dissTOM) <- yellow

# now choose the top connected genes for the genes of interest+ the top hub
yellow_genes_2_explore <- c('STARD8','POLN')
# finding the nearest neighbors of a gene is to look up the top items in sorted columns of adjacency
# matrix

yellow_adjacency <- adjacency[,yellow_genes_2_explore]
# stard8-first_gene
stard8 <- yellow_adjacency[order(yellow_adjacency$STARD8,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polb
poln <- yellow_adjacency[order(yellow_adjacency$POLN,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames

final_yellow <- c(stard8,poln)

# now to export module to cytoscape
modTOM = 1 - dissTOM
new_modTOM <- modTOM[final_yellow,final_yellow]
dimnames(new_modTOM) = list(final_yellow,final_yellow)
efname = "yellow_edges.txt"
nfname = "yellow_nodes.txt"

exportNetworkToCytoscape(new_modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = final_yellow,
)

names(yellow_nodes) <- c("nodeNames","altName","nodeAttr")
nodeAttr <- rep("yellow",22)
nodeAttr[1] <- "HUB"
nodeAttr[12] <- "TLS"
yellow_nodes$nodeAttr <- nodeAttr



