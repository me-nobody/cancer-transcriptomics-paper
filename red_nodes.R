
red_attr <- rep("red",22)
red_attr[12] <- "TLS"
red_attr[1] <- "HUB"
red_nodes$nodeAttr.nodesPresent... <- red_attr
write.csv(red_nodes,"red_nodes.txt",row.names = F,col.names = F)
