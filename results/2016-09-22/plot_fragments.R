L03  <- as.matrix(read.table('Lengths2.txt'))[,1]
L04  <- as.matrix(read.table('Lengths3.txt'))[,1]
L05  <- as.matrix(read.table('Lengths4.txt'))[,1]
L06  <- as.matrix(read.table('Lengths5.txt'))[,1]
L07  <- as.matrix(read.table('Lengths6.txt'))[,1]
L08  <- as.matrix(read.table('Lengths7.txt'))[,1]
L09  <- as.matrix(read.table('Lengths8.txt'))[,1]
L10  <- as.matrix(read.table('Lengths9.txt'))[,1]
L11  <- as.matrix(read.table('Lengths10.txt'))[,1]
L12  <- as.matrix(read.table('Lengths11.txt'))[,1]
L13  <- as.matrix(read.table('Lengths12.txt'))[,1]
L14  <- as.matrix(read.table('Lengths13.txt'))[,1]
L15  <- as.matrix(read.table('Lengths14.txt'))[,1]
L16  <- as.matrix(read.table('Lengths15.txt'))[,1]
L17  <- as.matrix(read.table('Lengths16.txt'))[,1]
L18  <- as.matrix(read.table('Lengths17.txt'))[,1]
L19  <- as.matrix(read.table('Lengths18.txt'))[,1]
L20  <- as.matrix(read.table('Lengths19.txt'))[,1]

Beta <- numeric(length=20)
Beta[3] <- mean(L03)
Beta[4] <- mean(L04)
Beta[5] <- mean(L05)
Beta[6] <- mean(L06)
Beta[7] <- mean(L07)
Beta[8] <- mean(L08)
Beta[9] <- mean(L09)
Beta[10] <- mean(L10)
Beta[11] <- mean(L11)
Beta[12] <- mean(L12)
Beta[13] <- mean(L13)
Beta[14] <- mean(L14)
Beta[15] <- mean(L15)
Beta[16] <- mean(L16)
Beta[17] <- mean(L17)
Beta[18] <- mean(L18)
Beta[19] <- mean(L19)
Beta[20] <- mean(L20)

K = 5.5
png(filename='QQ.png', width=6000, height=3000)
par(mfrow=c(3,6), mex=4, cex.main=K, cex.label=K, cex.axis=K, cex.lab=K)
qqplot(rexp(length(L03), rate=1/Beta[3]), L03, main=paste('Generation', 3), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L04), rate=1/Beta[4]), L04, main=paste('Generation', 4), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L05), rate=1/Beta[5]), L05, main=paste('Generation', 5), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L06), rate=1/Beta[6]), L06, main=paste('Generation', 6), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L07), rate=1/Beta[7]), L07, main=paste('Generation', 7), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L08), rate=1/Beta[8]), L08, main=paste('Generation', 8), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L09), rate=1/Beta[9]), L09, main=paste('Generation', 9), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L10), rate=1/Beta[10]), L10, main=paste('Generation', 10), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L11), rate=1/Beta[11]), L11, main=paste('Generation', 11), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L12), rate=1/Beta[12]), L12, main=paste('Generation', 12), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L13), rate=1/Beta[13]), L13, main=paste('Generation', 13), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L14), rate=1/Beta[14]), L14, main=paste('Generation', 14), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L15), rate=1/Beta[15]), L15, main=paste('Generation', 15), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L16), rate=1/Beta[16]), L16, main=paste('Generation', 16), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L17), rate=1/Beta[17]), L17, main=paste('Generation', 17), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L18), rate=1/Beta[18]), L18, main=paste('Generation', 18), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L19), rate=1/Beta[19]), L19, main=paste('Generation', 19), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(L20), rate=1/Beta[20]), L20, main=paste('Generation', 20), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
dev.off()

png(filename='FragmentLength.png', width=1000, height=1000)
par(mex=2, cex=2, cex.main=2, cex.axis=2, cex.lab=2)
plot(c(3:20), Beta[3:20]/100, xlab="Generations", ylab="Mean heterogenic fragment length (cM)")
dev.off()
