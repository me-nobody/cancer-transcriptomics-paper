library(tidyverse)
library(WGCNA)
tls_pol <- c("POLH","POLI","POLK","POLL","POLM","POLQ","POLZ","REV1","POLN")
repl_pol <- c("POLA1","POLB","POLG1","POLD","POLDE")
repair <- c("MLH1","MSH2","MSH6","BRCA2","RAD51","EXO1")
# brown,cyan,red,pink had relatively higher correlation with cetuximab
tls_pol%in%red
tls_pol%in%brown
tls_pol%in%pink
tls_pol%in%cyan
# above search showed no results and hence check to see if tls_pol
# is in the common_gene_subset
tls_pol%in%common_gene_subset
repair%in%common_gene_subset
repl_pol%in%common_gene_subset
# now search for modules where these are present
tls_pol_modules <- gene_color_df[gene_color_df$genes%in%tls_pol,]
repl_pol_modules <- gene_color_df[gene_color_df$genes%in%repl_pol,]
repair_pol_modules <- gene_color_df[gene_color_df$genes%in%repair,]
# the most predominant network is blue
# hence we shall create a network for blue for export to cytoscape

# first find genes in blue module
blue <- gene_color_df[gene_color_df$moduleColors=="blue",][,"genes"]
# subset combined_df by blue
blue_df <- corrected_df[,blue]
# softpower same as searching the original TOM for all genes
# apply soft-thresholding
softPower = 7
adjacency = adjacency(blue_df, power = softPower)

# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# now to export module to cytoscape
modTOM = 1 - dissTOM
dimnames(modTOM) = list(blue, blue)
efname = "blue_edges.txt"
nfname = "blue_nodes.txt"

exportNetworkToCytoscape(modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = blue,
                         )
# now we have to calculate module membership for all the genes in the blue
# module which is simply a  correlation of genes vs eigengene
cor_blue_genes<- cor(blue_df,MEs$MEblue)
blue_genes <- names(blue_df)
cor_blue_genes <- as.data.frame(cbind(blue_genes,cor_blue_genes))
names(cor_blue_genes) <- c("genes","mod.membership")
a <- cor_blue_genes[cor_blue_genes$genes%in%repl_pol,]
b <- cor_blue_genes[cor_blue_genes$genes%in%tls_pol,]
c <- cor_blue_genes[cor_blue_genes$genes%in%repair,]
blue_module_membership <- rbind(a,b,c)

# choose top hubs in each module in wgcna

hubs <- chooseTopHubInEachModule(corrected_df,moduleColors, omitColors = "grey")

# now choose the top connected genes for the genes of interest+ the top hub
blue_genes_2_explore <- c("RRM1",blue_module_membership$genes)

# finding the nearest neighbors of a gene is to look up the top items in sorted columns of adjacency
# matrix

blue_adjacency <- adjacency[,blue_genes_2_explore]

# rrmi-first_gene
rrmi <- blue_adjacency[order(blue_adjacency$RRM1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polb
polb <- blue_adjacency[order(blue_adjacency$POLB,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# pola1
pola1 <- blue_adjacency[order(blue_adjacency$POLA1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# poli
poli <- blue_adjacency[order(blue_adjacency$POLI,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polq
polq <- blue_adjacency[order(blue_adjacency$POLQ,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# polk
polk <- blue_adjacency[order(blue_adjacency$POLK,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# exo1
exo1 <- blue_adjacency[order(blue_adjacency$EXO1,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# brca2
brca2 <- blue_adjacency[order(blue_adjacency$BRCA2,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# rad51
rad51 <- blue_adjacency[order(blue_adjacency$RAD51,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# msh2
msh2 <- blue_adjacency[order(blue_adjacency$MSH2,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames
# msh6
msh6 <- blue_adjacency[order(blue_adjacency$MSH6,decreasing = T),1,drop=F]%>%head(n=11)%>%rownames

final_blue_genes <- c(rrmi,polb,pola1,poli,polq,polk,exo1,brca2,rad51,msh2,msh6)

# now to export module to cytoscape
modTOM = 1 - dissTOM
new_modTOM <- modTOM[final_blue_genes,final_blue_genes]
dimnames(new_modTOM) = list(final_blue_genes,final_blue_genes)
efname = "blue_edges.txt"
nfname = "blue_nodes.txt"

exportNetworkToCytoscape(new_modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.0, 
                         nodeNames = final_blue_genes,
)
