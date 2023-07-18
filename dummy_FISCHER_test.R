d <- data.frame(gene.not.interest=c(2613, 15310), gene.in.interest=c(28, 29))
row.names(d) <- c("In_category", "not_in_category")
fisher.test(d, alternative = "greater")
fisher.test(d)
library(enrichplot)
library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)