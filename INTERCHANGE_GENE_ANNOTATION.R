library(ggplot2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
# this library contains the relational db of gene ids and symbols
library(org.Hs.eg.db)
# use this argument to know the various key types for interconversion
keytypes(org.Hs.eg.db)
# the id list has to be a character vector with no special characters
# unquote removes quotes
blue_module <- noquote(blue_module)
# gsub will remove all special characters
blue_module <- lapply(blue_module, function(x) gsub('[^[:alnum:] ]','',x))
# this will convert the list to a character vector, else org.Hs.eg.db won't work
blue_module <- unlist(blue_module)
# to use enrichGO we need the entrez ids. we will use the select function
entrez_ids <- select(org.Hs.eg.db,keys=blue_module,keytype= "SYMBOL",columns = c("ENTREZID"))

edo <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
head(edo)
barplot(edo, showCategory=20) + theme(axis.text = element_text(size = 5))
