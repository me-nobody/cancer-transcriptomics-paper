# we have dataframes (SraRunTable) with Run information and Sample.Name. we have character array with Run names. we have
# to subset the dataframe based on the character array to get equivalent Sample.Names. then we have to read the text files
# corresponding to sample names and store the data in a dictionary with file names as keys and expression data as values.
library(tidyverse)
library(Dict)
library(stringr)
# for analysis select a subset of 50 runs for each dataset
GSE156451 <- sample(GSE156451_SraRunTable$Run,50)
GSE185055 <- sample(GSE185055_SraRunTable$Run,50)

# save the Run information as csv files

# filter sample_names corresponding to runs
GSE156451_Run_sample_match <- GSE156451_SraRunTable %>% filter(Run %in% GSE156451) 
GSE185055_Run_sample_match <- GSE185055_SraRunTable %>% filter(Run %in% GSE185055) 


# get corresponding matches with runs
gse156451_run_smpl_vec <- GSE156451_Run_sample_match$Sample.Name
gse185055_run_smpl_vec <- GSE185055_Run_sample_match$Sample.Name
# get the list of file names and change the variable file_list as required
file_list <- list.files(path="./datasets/GSE156451_RAW/")
file_list <- list.files(path="./datasets/GSE185055_RAW/")
# get shortened file names for file_list
file_list_truncated<- file_list %>% str_sub(1,10)
# get indexes in 2nd vector matching 1st vector
match_files<- match(gse156451_run_smpl_vec,file_list_truncated)
match_files<- match(gse185055_run_smpl_vec,file_list_truncated)
# get the desired files
desired_files <- file_list[match_files]
# create an empty dictionary. that init_key is needed else error is thrown
dataframe_dict <- Dict$new(init_keys = NULL)
folder_path <- "./datasets/GSE156451_RAW/"
folder_path <- "./datasets/GSE185055_RAW/"
for(file in desired_files){
      # collect file name as dictionary key
      entry_name <- str_sub(file,1,13)
      # set path to folder
      path <- paste(folder_path,file,sep="")
      # read the file as as a tibble and transform
      df <- read_tsv(path)%>%t
      colnames(df) <- df[1,]
      #df[1,] <- NULL
      dataframe_dict[entry_name] <- as_tibble(df)
      print(dim(dataframe_dict[entry_name]))
}
# this NULL key has to be removed now
dataframe_dict$remove("init_keys")
# once the NULL key is removed further operations can be done on the dataframe
for(key in dataframe_dict$keys){
  df <- dataframe_dict$get(key)
  df <- df[2,]
  row.names(df) <- key
  dataframe_dict[key] <- df
}
# create an emppty tibble with size and colnames as expression matrix
# df is from the earlier df created in the workspace
combined_df <- data.frame(matrix(nrow=0,ncol=length(df)))
names(combined_df) <- names(df)

for(key in dataframe_dict$keys){
  df <- dataframe_dict$get(key)
  combined_df <- rbind(combined_df,df)
}
# expression dataset formed
GSE156451_df <- combined_df
GSE185055_df <- combined_df
