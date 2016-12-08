lengths <- as.matrix(read.table('lengths.txt', header=TRUE))

png(filename='merged_length.png', width=1000, height=1000)
par(mfrow=c(2,2))
hist(lengths[,1], breaks=15, xlab="Length (bp)", main="Sample 1")
hist(lengths[,2], breaks=15, xlab="Length (bp)", main="Sample 3")
hist(lengths[,3], breaks=15, xlab="Length (bp)", main="Sample 5")
hist(lengths[,4], breaks=15, xlab="Length (bp)", main="Sample 6")
dev.off()
