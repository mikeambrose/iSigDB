args = commandArgs(TRUE)
sigFile = args[1]
samFile = args[2]
outFile = args[3]

sigs = read.table(sigFile,header=TRUE,sep="\t",row.names=1,check.names=FALSE)
datasets = read.table(samFile,header=TRUE,sep="\t",row.names=1,check.names=FALSE)

library(DeconRNASeq)

x = DeconRNASeq(datasets,sigs,fig=FALSE)
outMatrix = x$out.all
colnames(outMatrix) = colnames(sigs)
rownames(outMatrix) = colnames(datasets)
write.table(outMatrix,file=outFile,sep="\t",quote=FALSE)
