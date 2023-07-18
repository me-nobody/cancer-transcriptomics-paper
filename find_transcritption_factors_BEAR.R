BiocManager::install('KBoost')
# import the library for predicting transcription factors
library(KBoost)
library(Dict)
# list of TFs in KBoost database
human_tfs <- Human_TFs
# the dataset is expression matrix(NOT dataframe) with genes as columns and samples as rows
exprDat <- as.matrix(corrected_df)
# gene vector
genes <- names(corrected_df)
# calculate the transcription network using pre-trained features and default parameters
grn = KBoost_human_symbol(exprDat,genes,pos_weight=0.6, neg_weight=0.4)
prob_matrix_grn <- grn$GRN
prior_grn <- grn$prior_weights
TF_list <- grn$TFs
grn_TFs <- grn$model
# transcription factors in the corrected_df expression dataset
TFs_in_corrected_df <- genes[TF_list]

# list of TFs in KBoost database
human_tfs <- Human_TFs
# check validity ofTFs_in_corrected_df
all(TFs_in_corrected_df%in%human_tfs)
# names of TFs and genes in grn_TFs 
names(grn_TFs) <- TFs_in_corrected_df
rownames(grn_TFs) <- names(corrected_df)
# select target size as  30 to screen for TFs
logic_vec <- apply(grn_TFs,2,function(x) sum(x)>=30)
# find the TFs with 30 or more targets
filtered_tfs <- TFs_in_corrected_df[logic_vec]  
# filter the predicted TF dataframe with the filtered TF list
grn_TFs_filtered <- grn_TFs[filtered_tfs]
# combine the list of key genes in the modules
key_genes<- c(final_blue_genes,final_greenyellow,final_red,final_turquoise,final_yellow)
# filter the TF dataframe by key_genes
grn_TFs_filtered<- grn_TFs_filtered[rownames(grn_TFs_filtered)%in%key_genes,]
# further filtering TFs which act on the key genes
TFs_final_filter <- names(grn_TFs_filtered)[apply(grn_TFs_filtered,2,sum)>0]
# subset grn_TFs_filtered based on above filter
grn_TFs_filtered <- grn_TFs_filtered[TFs_final_filter]
tfs_targets <- Dict$new(
   dummy = '',
  .overwrite = TRUE
)

for (col in seq_along(grn_TFs_filtered)){
  filter_query <- which(grn_TFs_filtered[col]==TRUE)
  targets <- row.names(grn_TFs_filtered)[filter_query]
  tfs_targets[names(grn_TFs_filtered[col])]=targets
  
  }

file_conn = file("TFs_n_targets.txt")
close(file_conn)
for (item in tfs_targets$keys){
  cat("**************",file = "TFs_n_targets.txt",append=TRUE,sep="\n")
  cat("          ",sep = "\n",append = TRUE)
  cat(paste("Transcripton Factor",item,sep = "\t"),file="TFs_n_targets.txt",append=TRUE,sep = "\n")
  cat("--------",file = "TFs_n_targets.txt",append=TRUE,sep = "\n")
  cat("          ",sep = "\n",append = TRUE)
  cat(dummy[item],file= "TFs_n_targets.txt",append=TRUE,sep="\n")}
  cat("          ",sep = "\n",append = TRUE)
file.show("TFs_n_targets.txt")

# function to calculate shortest distance between nodes
dist_0 <- net_dist_bin(GRN,TFs,thr)





