bamsize <- read.table('bamsize.txt', header=TRUE)
attach(bamsize)

png(filename='bamsize.png')
plot(FastqSize, BamSize, xlab="Fastq size (bytes)", ylab="Bam size (bytes)")
dev.off()

detach(bamsize)
