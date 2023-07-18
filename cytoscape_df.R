library(WGCNA)

# create the new dataframe for the chosen genes from the chosene modules
blue_cytoscape <- final_blue_genes
red_cytoscape <- final_red
greenyellow_cytoscape <- final_greenyellow
turquoise_ctyoscape <- final_turquoise
yellow_cytoscape <- final_yellow
# combined  vector of the desired genes for adjacency matrix
key_genes <- c(blue_cytoscape,red_cytoscape,greenyellow_cytoscape,turquoise_ctyoscape,yellow_cytoscape)

# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency_cytoscape = adjacency(corrected_df, power = softPower)
adjacency <- adjacency_cytoscape[key_genes,key_genes]
# create topological overlap matrix
TOM_cytoscape = TOMsimilarity(adjacency)
dissTOM_cytoscape = 1-TOM_cytoscape

colnames(TOM_cytoscape) <- key_genes
rownames(TOM_cytoscape) <- key_genes
colnames(dissTOM_cytoscape) <- key_genes
rownames(dissTOM_cytoscape) <- key_genes
# now to export module to cytoscape
# for cytoscape you basically need the TOM file

exportNetworkToCytoscape(TOM_cytoscape,
                         edgeFile = "cytoscape_edges.txt",
                         nodeFile = "cytoscape_nodes.txt",
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = key_genes,
)
