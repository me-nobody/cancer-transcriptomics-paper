library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

blue_module <- noquote(blue_module)
blue_module <- lapply(blue_module, function(x) gsub('[^[:alnum:] ]','',x))
blue_module <- unlist(blue_module)
edo <- enrichGO(gene          = blue_module,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
barplot(edo, showCategory=20) 