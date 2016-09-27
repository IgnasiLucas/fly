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

E03  <- as.matrix(read.table('equalRates/Lengths2.txt'))[,1]
E04  <- as.matrix(read.table('equalRates/Lengths3.txt'))[,1]
E05  <- as.matrix(read.table('equalRates/Lengths4.txt'))[,1]
E06  <- as.matrix(read.table('equalRates/Lengths5.txt'))[,1]
E07  <- as.matrix(read.table('equalRates/Lengths6.txt'))[,1]
E08  <- as.matrix(read.table('equalRates/Lengths7.txt'))[,1]
E09  <- as.matrix(read.table('equalRates/Lengths8.txt'))[,1]
E10  <- as.matrix(read.table('equalRates/Lengths9.txt'))[,1]
E11  <- as.matrix(read.table('equalRates/Lengths10.txt'))[,1]
E12  <- as.matrix(read.table('equalRates/Lengths11.txt'))[,1]
E13  <- as.matrix(read.table('equalRates/Lengths12.txt'))[,1]
E14  <- as.matrix(read.table('equalRates/Lengths13.txt'))[,1]
E15  <- as.matrix(read.table('equalRates/Lengths14.txt'))[,1]
E16  <- as.matrix(read.table('equalRates/Lengths15.txt'))[,1]
E17  <- as.matrix(read.table('equalRates/Lengths16.txt'))[,1]
E18  <- as.matrix(read.table('equalRates/Lengths17.txt'))[,1]
E19  <- as.matrix(read.table('equalRates/Lengths18.txt'))[,1]
E20  <- as.matrix(read.table('equalRates/Lengths19.txt'))[,1]

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

Veta <- numeric(length=20)
Veta[3] <- mean(E03)
Veta[4] <- mean(E04)
Veta[5] <- mean(E05)
Veta[6] <- mean(E06)
Veta[7] <- mean(E07)
Veta[8] <- mean(E08)
Veta[9] <- mean(E09)
Veta[10] <- mean(E10)
Veta[11] <- mean(E11)
Veta[12] <- mean(E12)
Veta[13] <- mean(E13)
Veta[14] <- mean(E14)
Veta[15] <- mean(E15)
Veta[16] <- mean(E16)
Veta[17] <- mean(E17)
Veta[18] <- mean(E18)
Veta[19] <- mean(E19)
Veta[20] <- mean(E20)

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

png(filename='equalRates/QQ.png', width=6000, height=3000)
par(mfrow=c(3,6), mex=4, cex.main=K, cex.label=K, cex.axis=K, cex.lab=K)
qqplot(rexp(length(E03), rate=1/Veta[3]), E03, main=paste('Generation', 3), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E04), rate=1/Veta[4]), E04, main=paste('Generation', 4), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E05), rate=1/Veta[5]), E05, main=paste('Generation', 5), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E06), rate=1/Veta[6]), E06, main=paste('Generation', 6), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E07), rate=1/Veta[7]), E07, main=paste('Generation', 7), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E08), rate=1/Veta[8]), E08, main=paste('Generation', 8), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E09), rate=1/Veta[9]), E09, main=paste('Generation', 9), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E10), rate=1/Veta[10]), E10, main=paste('Generation', 10), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E11), rate=1/Veta[11]), E11, main=paste('Generation', 11), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E12), rate=1/Veta[12]), E12, main=paste('Generation', 12), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E13), rate=1/Veta[13]), E13, main=paste('Generation', 13), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E14), rate=1/Veta[14]), E14, main=paste('Generation', 14), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E15), rate=1/Veta[15]), E15, main=paste('Generation', 15), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E16), rate=1/Veta[16]), E16, main=paste('Generation', 16), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E17), rate=1/Veta[17]), E17, main=paste('Generation', 17), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E18), rate=1/Veta[18]), E18, main=paste('Generation', 18), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E19), rate=1/Veta[19]), E19, main=paste('Generation', 19), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
qqplot(rexp(length(E20), rate=1/Veta[20]), E20, main=paste('Generation', 20), cex=K, xlab="Exponential quantiles", ylab="Observed quantiles")
abline(0, 1, col="red", lwd=K)
dev.off()


png(filename='FragmentLength.png', width=1000, height=1000)
par(mex=2, cex=2, cex.main=2, cex.axis=2, cex.lab=2)
plot(c(3,20),c(0, max(Beta,Veta)/100), type="n", xlab="Generations", ylab="Mean heterogenic length (cM)")
lines(3:20, Beta[3:20]/100, col="red", lwd=3)
lines(3:20, Veta[3:20]/100, col="blue", lwd=3)
dev.off()
