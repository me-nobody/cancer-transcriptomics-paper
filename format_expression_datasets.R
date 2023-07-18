library(dplyr)
library(biomaRt)
GSE196576_df <- GSE196576_CALGB_80405_RNAseq_Processed_data
GSE196576_df <- GSE196576_df[,-c(1,3:8)]
GSE196576_df <- as.data.frame(t(GSE196576_df))
names(GSE196576_df) <- GSE196576_df[1,]
GSE196576_df <- GSE196576_df[-1,]
# Accession has NA values and hence will not be accepted as rownames, so we will keep the original names
row.names(GSE196576_df) <- GSE196576_df$Accession
which(is.na(GSE196576_df$Accession))
# we will proceed with custom values provided by the PI as Accession has
# NA values and we need to proceed with the Cetuximab subset
GSE196576_Cetuximab <- GSE196576_Cetuximab%>%left_join(GSE196576_accession_with_sample_names,by="Accession")
# title column of GSE196576_Cetuximab
title_cetuximab <- GSE196576_Cetuximab$Title
# get the subset of GSE196576 corresponding to Cetuximab df
GSE196576_cetuximab <- GSE196576_df[rownames(GSE196576_df)%in%title_cetuximab,]
# get a sample of 50 from the subset
GSE196576_cetuximab <- GSE196576_cetuximab[sample(rownames(GSE196576_cetuximab),50),]
#------------------------------------------------------------
# GSE185055 has ensemble ids. we have to convert it to gene names
# we have downloaded the GRch19 ID gene table from ENSEMBL and not using biomart for now
names_gse185055_df <- as.data.frame(names(GSE185055_df))
names(names_gse185055_df) <- "Gene.stable.ID"
#human_GRch37_ensembl_gene_conversion is the file downloaded with the gene name ENSEMBLE matching 
names_gse185055_ided <- names_gse185055_df%>%left_join(human_GRch37_ensembl_gene_conversion,by="Gene.stable.ID")
# apply the gene names to the GSE185055_df dataframe
names(GSE185055_df) <- names_gse185055_ided$Gene.name
#--------------------------------------------------------------
# get a common list of genes GSE196576/GSE156451/GSE185055
#--e.g for multiple intersections intersect(intersect(a,b),c)
# common genes 
a <- names(GSE156451_df)
b <- names(GSE196576_cetuximab)
c <- names(GSE185055_df)
common_genes_cetuximab_tumor_celline <- intersect(intersect(a,b),c)
# reset GSE196576 and GSE156451 with common set of genes
GSE196576_expr_cetuximab <- GSE196576_cetuximab[,common_genes_cetuximab_tumor_celline]
GSE156451_expr_tumor <- GSE156451_df[,common_genes_cetuximab_tumor_celline]
rownames(GSE156451_expr_tumor) <- rownames(GSE156451_df)
# apply the same subset for GSE185055 as the other dataframes
GSE185055_expr_celline <- GSE185055_df[,common_genes_cetuximab_tumor_celline]
rownames(GSE185055_expr_celline) <- row.names(GSE185055_df)
GSE156451_expr_tumor <- as.data.frame(GSE156451_expr_tumor)

GSE185055_expr_celline <- as.data.frame(GSE185055_expr_celline)

#---------------------------------------------------
#-----------------discarded steps-------------------
# creating a file with GSM as sample ids
GSE196576_accession_with_sample_names <- GSE196576_accession_with_sample_names[,c("Accession","Title","SRA.Accession")]
GSE196576_df <- cbind(row.names(GSE196576_df),GSE196576_df)
names(GSE196576_df) <- c("Title",names(GSE196576_df)[-1])

# remove all columns with NA else left join doesn't work
GSE196576_df <- GSE196576_df[,!(is.na(names(GSE196576_df)))]
# add accession to the dataframe
GSE196576_df <- GSE196576_df%>%left_join(GSE196576_accession_with_sample_names,by="Title")
