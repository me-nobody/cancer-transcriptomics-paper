library(WGCNA)
library(dplyr)

# transpose the datasets- samples as rows and genes as columns
combat_edate_par <- t(combat_edate_par)

# search soft-thresholding powers
powers = 2:20
sft = pickSoftThreshold(combat_edate_par, powerVector = powers, verbose = 5)

# plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2', type = 'n',
     main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = 'red');
abline(h = 0.90, col = 'red')
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = 'red')


net = blockwiseModules(combat_edate_par, power = 7,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "merged_1st_network",
                          verbose = 1)

table(net$colors)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
summary(as.factor(moduleColors))
MEs = net$MEs
corrplot(cor(MEs))
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "first_network-auto.RData")
# Define numbers of genes and samples
nGenes = ncol(combat_edate_par);
nSamples = nrow(combat_edate_par);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(combat_edate_par, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
corrplot(cor(MEs0))
