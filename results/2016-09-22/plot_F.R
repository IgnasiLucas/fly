library(RColorBrewer)

ibd <- read.table('F.txt', header=TRUE)
attach(ibd)

# Theoretical expectations of inbreeding coefficient under full-sib mating are
# taken from Lynch and Walsh (1998), page 259, equation 10.5b.

TheoreticalF <- numeric(length = 20)
TheoreticalF[1] = 0
TheoreticalF[2] = 0.25
for (i in 3:20) {
   TheoreticalF[i] <- TheoreticalF[i-1]/2 + TheoreticalF[i-2]/4 + 1/4
}

AverageF <- aggregate(F1, by=list(Gen), mean)

png(filename="F.png", width=1000, height=1000)
par(mex=2, cex.axis=2, cex.lab=2)
boxplot(F1 ~ Gen, xlab="Generations", ylab="Inbreeding coefficient")
lines(c(1:20), AverageF[,2], col="red", lwd = 2)
lines(c(1:20), TheoreticalF, col="blue", lwd = 2)
dev.off()


png(filename="NumberTracts.png", width=1000, height=1000)
par(mex=2, cex.axis=2, cex.lab=2, cex.main=2)
barplot(table(N.Tracts.1, Gen), col=brewer.pal(7,'YlOrRd'), xlab="Generations", ylab="Frequency", main="Number of heterogenic tracts")
dev.off()
