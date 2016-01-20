if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

options(scipen=10)
args<-commandArgs(TRUE)
file=args[1]
file_label=gsub(".fasta_merged.txt_toR.txt","",file)
par(ps = 22, cex = 0.7, cex.main = 0.4)

#heatmap by genes
data = read.table(file, sep=" ", header=TRUE, stringsAsFactors=F)
data_transformed=cbind(data[, 1], round(
        data[, seq(from = 3, to = ncol(data), by = 3)] *100 /
        data[, seq(from = 2, to = ncol(data), by = 3)], digits=2))
data_transformed[data_transformed==Inf] <- 0
row.names(data_transformed) <- data_transformed[,1]
labels=data_transformed[,1]
data_transformed=data_transformed[,-1]
data_before=data_transformed 
#recalculate over represented genes (coverage greater than 100)
over_covered <- data_transformed > 100

#USE THIS WHEN INTERESTED IN PENALIZING OVERCOVERED REGIONS
#data_transformed[over_covered] <- data_transformed[over_covered] - 2 * (data_transformed[over_covered] - 100)

#USE THIS WHEN OVERCOVERED REGIONS SHOULD BE SET TO 100 %
data_transformed[over_covered] <- 100

#complete data
data_complete=cbind(labels,data_transformed)
data_complete <- data_complete[,2:ncol(data_complete)]
data_matrix_genes <- data.matrix(data_complete)
data_matrix_genes_original_values <- data.matrix(data_before)

mycol<-colorRampPalette(c("#fee6ce","#fdae6b","#e6550d"), bias=1)

outputfile=paste('figures/','heatmapByGenes_',file_label,'.pdf',sep = "");
pdf(outputfile, height=12, width=12)
data_heatmap <- heatmap.2(data_matrix_genes, Rowv=NA, Colv=NA, col=mycol(4), scale="none", 
                          key=FALSE, trace="none", #tracecol="#054662",
                          dendrogram="none", lhei=c(1,12), cexRow=1, cexCol=1,
                          margins = c(12, 8), offsetCol = 0.01, main = file_label, 
                          sepwidth=c(0.000001,0.000005),
                          colsep=0:ncol(data_matrix_genes),
                          rowsep=0:nrow(data_matrix_genes),
                          breaks=c(0,25,50,75,100),
                          )

dev.off()