library(KBoost)
exprMat <- t(blue_df)
set.seed(123)
out <- KBoost_human_symbol(exprMat)
