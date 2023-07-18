# importing required libraries
library(tidyverse) # bullshit library for melting
library(reshape2)
library(ggplot2)
# calculate the z-score
z_dist_chosen_tf_targets <- apply(dist_chosen_tf_targets,2,function(x) (0.3-mean(x))/sd(x))
# creating the modified dataframe
# we loose the information on the targets
tf_dist_long <- dist_chosen_tf_targets
tf_dist_long['targets'] <- row.names(tf_dist_long)
tf_dist_long%>%relocate(targets)
tf_dist_long <- melt(tf_dist_long,variable.name='TFs',value.name='dist',id.vars='targets')

# creating a plot
p <- ggplot(tf_dist_long) + geom_boxplot(aes(x=TFs, y=dist, color=TFs))

# printing the plot
plot(p)
