library(Hmisc)
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(correlation_matrix, p_matrix) {
  ut <- upper.tri(correlation_matrix)
  data.frame(
    row = rownames(correlation_matrix)[row(correlation_matrix)[ut]],
    column = rownames(correlation_matrix)[col(correlation_matrix)[ut]],
    cor  =(correlation_matrix)[ut],
    p = p_matrix[ut]
  )
}

MEcorr<-rcorr(as.matrix(MEs))
MEcorr_df <- flattenCorrMatrix(MEcorr$r, MEcorr$P)
