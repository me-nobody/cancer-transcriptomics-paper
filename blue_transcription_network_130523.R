library(tidyverse)
blue_nodes_updated <- blue_nodes%>%left_join(tf_target_nodes,by=c('nodeName'='name'))
blue_nodes_updated <- blue_nodes_updated%>%left_join(tf_target_edge,by=c('nodeName'='targets'))
write.table(blue_nodes_updated,"blue_nodes_tf_updated.txt",
            row.names = F,quote = F,sep="\t")
