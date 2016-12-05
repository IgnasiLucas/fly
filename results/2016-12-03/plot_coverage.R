pooled <- read.table('coverage_pooled.txt', col.names=c("cov", "freq", "reads"), header=TRUE)
i1b1   <- read.table('coverage_i1b1.txt', col.names=c("cov", "freq", "reads"), header = TRUE)
i3b3   <- read.table('coverage_i3b3.txt', col.names=c("cov", "freq", "reads"), header = TRUE)
i5b5   <- read.table('coverage_i5b5.txt', col.names=c("cov", "freq", "reads"), header = TRUE)
i6b6   <- read.table('coverage_i6b6.txt', col.names=c("cov", "freq", "reads"), header = TRUE)


png(filename="coverage.png", height=500, width=1000)
par(mfrow=c(1,2))
plot(c(0,500), c(0,200), type='n', xlab='Coverage', ylab='Number of loci')
lines(i1b1$cov, i1b1$freq, col='red')
lines(i3b3$cov, i3b3$freq, col='orange')
lines(i5b5$cov, i5b5$freq, col='green')
lines(i6b6$cov, i6b6$freq, col='blue')

plot(c(0,1000), c(0,1), type='n', xlab='Coverage', ylab='Cum. Prop. of reads')
lines(i1b1$cov, i1b1$reads/i1b1$reads[length(i1b1$reads)], col='red')
lines(i3b3$cov, i3b3$reads/i3b3$reads[length(i3b3$reads)], col='orange')
lines(i5b5$cov, i5b5$reads/i5b5$reads[length(i5b5$reads)], col='green')
lines(i6b6$cov, i6b6$reads/i6b6$reads[length(i6b6$reads)], col='blue')

dev.off()
