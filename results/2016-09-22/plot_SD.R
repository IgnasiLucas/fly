ibdm0 <- read.table('F.txt', header=TRUE)
ibdeq <- read.table('equalRates/F.txt', header=TRUE)

StdDevFm0 <- aggregate(ibdm0$F1, by=list(ibdm0$Gen), sd)
StdDevFeq <- aggregate(ibdeq$F1, by=list(ibdeq$Gen), sd)

png(filename="SD.png", width=1000, height=1000)
par(mex=2, cex.axis=2, cex.lab=2)
plot(StdDevFm0[,1], StdDevFm0[,2], type="l", lwd=2, xlab="Generations", ylab="Standard Deviation of F", col="red")
lines(StdDevFeq[,1], StdDevFeq[,2], col="blue", lwd=2)
dev.off()
