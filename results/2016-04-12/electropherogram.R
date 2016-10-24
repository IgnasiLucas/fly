bp <- as.matrix(read.table('RCATGY/bp.txt'))
low        <- 20
while (is.na(match(low, bp[,1]))) {
   low <- low + 1
}
high       <- 10500
while (is.na(match(high, bp[,1]))){
   high <- high + 1
}
index.low  <- match(low, bp[,1])
index.high <- match(high, bp[,1])
value.low  <- min(bp[index.low:index.high, 2])
value.high <- max(bp[index.low:index.high, 2])
png(filename = 'electropherogram.png', width = 1000, height = 1000)
par(mex = 2)
plot(c(low, high), c(value.low, value.high), type="n", log="x", 
   cex.main = 2, cex.axis = 2, cex.lab = 2, main = "Expected electropherogram",
   xlab = "Fragment length (bp)", ylab = "Accumulated number of bp")
lines(bp[index.low:index.high, 1], bp[index.low:index.high, 2], lwd = 2)
dev.off()
