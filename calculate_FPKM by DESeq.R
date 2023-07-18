library(dplyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
par(mfrow=c(3,1))
tumor <- sapply(GSE156451,sum)
cellline <- sapply(GSE185055,sum)
cetuximab <- sapply(GSE196576,sum)
hist(sapply(tumor,log2),main="raw counts tumor",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
hist(sapply(cellline,log2),main="raw counts cell line",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
hist(sapply(cetuximab,log2),main="raw counts cetuximab",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
dev.off()

# we have again lost the rownames, so proceedign w/o them
rownames(GSE156451) <- paste("tumor",1:50,sep = "_")
rownames(GSE185055) <- paste("cellline",1:50,sep = "_")
rownames(GSE196576) <- paste("cetuximab",1:50,sep = "_")

# we can't combine right away as GSE156451 is FPKM while GSE185055,
# and GSE196576 are count data
