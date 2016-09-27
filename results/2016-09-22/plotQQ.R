gen <- as.numeric(commandArgs(trailingOnly=TRUE))[1]
L   <- as.matrix(read.table(paste0('Lengths', gen, '.txt')))[,1]
Beta <- mean(L)
expo <- rexp(length(L), rate=1/Beta)
png(filename=paste0('Lengths',gen,'.png'), width=1000, height=1000)
par(mex=2, cex.main=2, cex.label=2, cex.axis=2, cex.lab=2)
qqplot(expo, L, main=paste('Generation', gen), cex=2)
abline(0, 1, col="red", lwd=2)
dev.off()
