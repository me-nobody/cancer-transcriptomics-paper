library(WGCNA)
library(tidyverse)
library(corrplot)
# transpose the datasets- samples as rows and genes as columns
corrected_df <- as.data.frame(t(corrected_df))

# search soft-thresholding powers
powers = 2:20
sft = pickSoftThreshold(corrected_df, powerVector = powers, verbose = 5)

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


net = blockwiseModules(corrected_df, power = 7,
                          networkType = "signed",minKMEtoStay = 0.5,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.5,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          saveTOMFileBase = "tls_pol_network_080523",
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
corrplot(cor(MEs),tl.cex = 0.7)
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "second_network-auto.RData")
# Define numbers of genes and samples
nGenes = ncol(corrected_df);
nSamples = nrow(corrected_df);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(corrected_df, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# remove the grey module as these genes could not be networked
# have to use full dplyr path as select function overlaps with annotationdbi
MEs <- dplyr::select(MEs,-MEgrey)
# plotting high resolution for poster
tiff("Eigengene_correlation.tiff", width = 6, height = 6, units = 'in', res = 300)
corrplot(cor(MEs),tl.cex = 0.9,tl.srt=30,
         diag=FALSE,order="hclust",type="lower",addCoef.col = 'black',
         number.cex=0.6)
dev.off()
par()$pin # plot parameters

# trait_module_correlation
# Define a categorical variable with 3 levels
x = rep(c("tumor", "cell_line", "cetuximab"), each = 50)
# Binarize it into pairwise indicators
out = binarizeCategoricalVariable(x,
                                  includePairwise = FALSE,
                                  includeLevelVsAll = TRUE);
# Print the variable and the indicators
out <- data.frame(out)
names(out) <- c("cell_line","cetuximab","tumor")
traits <- out%>%relocate(tumor,.before = cell_line)
rownames(traits) <- rownames(corrected_df)
# traits dataframe is ready
# sample traits to have same rows as MEs, which is sample names

# actual correlation calculation
moduleTraitCor = cor(MEs, traits, method = 'kendall');
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(traits));

# plot the results
tiff("CRC_kendall_correlation_module_trait.tiff", width = 6, height = 6, units = 'in', res = 300)
sizeGrWindow(10,10)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 12, 4, 12));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.axis = 0.4,
               cex.lab =0.6,
               zlim = c(-1,1),
               main = paste("CRC transcriptome trait-module relation"))

dev.off()

# poor correlation with traits-attempting glm
# glm for tumor samples
glm_tumor_pval <- apply(MEs,2,function(x)  summary(glm(traits$tumor~x,
                                family = 'binomial'(link = 'logit')))$coefficients[8])

# glm for cell_lines samples
glm_cell_line_pval <- apply(MEs,2,function(x)  summary(glm(traits$cell_line~x,
                                            family = 'binomial'(link = 'logit')))$coefficients[8])

# glm for cetuximab samples
glm_cetuximab_pval <- apply(MEs,2,function(x)  summary(glm(traits$cetuximab~x,
                                              family = 'binomial'(link = 'logit')))$coefficients[8])

glm_res_df <- as.data.frame(rbind(glm_tumor_pval,glm_cell_line_pval,glm_cetuximab_pval))
# kendall correlation better than glm fit
# script to extract genes but we will use a different approach here.
grab_genes <- function(module_list, data) {
  mod_list <- character()
  for(color in module_list){
    genes <- colnames(data)[moduleColors == color]
    mod_list <- c(mod_list,genes)
    return (mod_list)
  }
}

gene_color_df <- as.data.frame(cbind(moduleLabels,moduleColors))
gene_color_df$moduleLabels <- as.integer(gene_color_df$moduleLabels)

# these modules were statistically significant by module trait relationship
pink <- rownames(gene_color_df[gene_color_df$moduleColors=="pink",])
red <- rownames(gene_color_df[gene_color_df$moduleColors=="red",])
cyan <- rownames(gene_color_df[gene_color_df$moduleColors=="cyan",])
brown <- rownames(gene_color_df[gene_color_df$moduleColors=="brown",])
write.csv(pink,"pink_module.csv",row.names = F,col.names = F)
write.csv(red,"red_module.csv",row.names = F,col.names = F)
write.csv(cyan,"cyan_module.csv",row.names = F,col.names = F)
write.csv(brown,"brown_module.csv",row.names = F,col.names = F)