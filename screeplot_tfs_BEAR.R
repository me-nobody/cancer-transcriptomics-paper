screeplot_tfs <- function(TFs_by_KBoost,TFs_in_expression_data) {
               num_targets <- 1:50
               num_tfs <- c()
                    for(i in num_targets){
                    logic_vec <- apply(grn_TFs,2,function(x) sum(x)>i)
                    filtered_tfs <- TFs_in_corrected_df[logic_vec]  
                    if(i==30){print(length(filtered_tfs))}
                    num_tfs <- c(num_tfs,length(filtered_tfs))
                    }
               par(oma=c(3,2,2,2))
               plot(num_tfs,num_targets,xlab="number of TFs",
                    ylab="number of targets",main="TF target prediction by KBoost",
                    cex=1.5,cex.lab=2,cex.axis=2,pch=16)
               abline(h=30,v=144,col='red',lwd=4)
                }
screeplot_tfs(grn_TFs,TFs_in_corrected_df)
par(oma=c(3,3,2,2))
