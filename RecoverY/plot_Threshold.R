h1k <- read.delim("dsk_output/head1k_kmers", header=FALSE)
jpeg("myplot.jpg")
plot (h1k$V1, log(h1k$V2), col="orange", ylab='Number of k-mers(logarithmic scale)', xlab='Number of times a distinct k-mer occurs')
dev.off()
