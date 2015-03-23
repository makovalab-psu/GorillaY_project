if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

options(scipen=10)
#setwd("/Users/alice/Desktop/projects/dotplot/heatmaps")
args<-commandArgs(TRUE)
file=args[1]
#file="P6-P8_merged.txt"

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
#over_covered[,2:ncol(data_transformed)] <- data_transformed[,2:ncol(data_transformed)] > 100
over_covered <- data_transformed > 100
#View(over_covered)


#USE THIS WHEN INTERESTED IN PENALIZING OVERCOVERED REGIONS
#data_transformed[over_covered] <- data_transformed[over_covered] - 2 * (data_transformed[over_covered] - 100)

#USE THIS WHEN OVERCOVERED REGIONS SHOULD BE SET TO 100 %
data_transformed[over_covered] <- 100

#View(data_transformed)

#complete data
data_complete=cbind(labels,data_transformed)
#View(data_complete)
#View(data_before)

data_complete <- data_complete[,2:ncol(data_complete)]
#data_complete$avg <- round(rowMeans(data_complete, na.rm = TRUE), digits=2)
#data_complete <- data_complete[order(data_complete$avg),]

#data_before <- data_before[,2:ncol(data_before)]
#data_before$avg <- round(rowMeans(data_before, na.rm = TRUE), digits=2)
#data_before <- data_before[order(data_complete$avg),]

#plot(data_complete$avg, data$length)

data_matrix_genes <- data.matrix(data_complete)
data_matrix_genes_original_values <- data.matrix(data_before)

mycol<-colorRampPalette(c("red","purple","darkblue"))(10)

outputfile=paste('figures/','heatmapByGenes_',file_label,'.pdf',sep = "");
pdf(outputfile, height=12, width=12)
data_heatmap <- heatmap.2(data_matrix_genes, Rowv=NA, Colv=NA, col=mycol, scale="none", 
                          key=FALSE, trace="none", tracecol="black",
                          dendrogram="none", lhei=c(1,12), cexRow=1, cexCol=1,
                          margins = c(12, 8), offsetCol = 0.01, main = file_label, 
                          sepwidth=c(0.00001,0.00001),
                          sepcolor="white",
                          colsep=1:ncol(data_matrix_genes),
                          rowsep=1:nrow(data_matrix_genes),
                          breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                          cellnote=data_complete,
                          notecol="white"
                          )

dev.off()