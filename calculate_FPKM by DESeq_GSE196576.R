library(dplyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# first transpose GSE196576
GSE196576 <- as.data.frame(t(GSE196576))
sampleNames <- names(GSE196576)
condition <- as.factor(rep("cetuximab",50))
colData <- as.data.frame(cbind(sampleNames,condition))
colData$sampleNames==names(GSE196576)
# create the deseq object
dds_gse196576 <- DESeqDataSetFromMatrix(countData = GSE196576,colData = colData, design = ~1 )
# here design=~1 has to be used as this is not a differential expression issue.
dds_gse196576 <- DESeq(dds_gse196576)
# filter the low count genes to remove noise
keep_gse196576 <- rowSums(counts(dds_gse196576))>=10
dds_76 <- dds_gse196576[keep_gse196576,]
# the dataframe gene_name_gene_lengths was obtained from ENSEMBL biomart.
# we have to derive the basepairs column to provide as input for FPKM
gene_lengths <- gene_name_gene_lengths
names(gene_lengths) <- c("genes","start","end")
gene_lengths$basepairs <- gene_lengths$end-gene_lengths$start
# get the list and order of genes in the dds object
dds_genes_gse196576 <- rownames(dds_76)
# use dplyr filter to get the genes common
gene_lengths_gse196576 <- gene_lengths%>%filter(genes%in%dds_genes_gse196576)
# there are duplicates and hence use dplyr distinct to get unique entries
gene_lengths_gse196576 <- gene_lengths_gse196576%>%distinct(genes,.keep_all = TRUE)
# now we have to reorder the gene names as per the dds object to get the gene lengths
# first create a dataframe from dds_genes_gse196576
dds_genes_gse196576 <- as.data.frame(dds_genes_gse196576)
names(dds_genes_gse196576) <- "genes"
# now use dplyr left_join to merge gene_lengths_gse196576 with this dataframe
dds_genes_gse196576 <- dds_genes_gse196576%>%left_join(gene_lengths_gse196576,by="genes")
# now from the above dataframe we  will get the genelengths to calculate FPKM
# we will assign these gene lengths to GSE196576
mcols(dds_76)$basepairs <- dds_genes_gse196576$basepairs
# now calculation of FPKM by DESeq2
fpkm_gse196576 <- as.data.frame(fpkm(dds_76))



keep_gse185055 <- rowSums(counts(dds_gse185055))>=10
dds_55 <- dds_gse185055[keep_gse185055,]

