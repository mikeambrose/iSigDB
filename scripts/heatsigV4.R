## heatSig.R
#==============================================================================
# Date Modified: 2015.04.24
#==============================================================================
args = commandArgs(TRUE)
srcFile = args[1]
outFile = args[2]
bcolumn = args[3]
rowClustMethod = args[4]
colClustMethod = args[5]
outTxtFile = args[6]

library(gplots);
rma = read.table(srcFile,header=TRUE,sep = "\t", row.names=1, comment.char="",check.names=FALSE)
numRows = nrow(rma)
numCols = ncol(rma)
rmamat = as.matrix(rma)

rowClust = (rowClustMethod != "none") && (numRows != 1)
colClust = (colClustMethod != "none") && (numCols != 1)

if (bcolumn=="matrix") {
  matMean <- mean(as.vector(rmamat))
  matSD <- sd(as.vector(rmamat))
  for (row in 1:nrow(rmamat)) {
    for (col in 1:ncol(rmamat)) {
      originalValue <- rmamat[row, col];
      rmamat[row, col] = (originalValue - matMean) / (matSD);
    }
  }
}

if (rowClust) {

  if (rowClustMethod == "pear_cor") {
    hr <- hclust(as.dist(1-cor(t(rmamat), method="pearson")), method="complete");
  } else if (rowClustMethod == "euclidean") {
    hr <- hclust(dist(rmamat, method="euclidean"), method="complete");
  }

  mycl <- cutree(hr, h=max(hr$height)/12);
  mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9);
  mycolhc <- mycolhc[as.vector(mycl)]
}
if (colClust) {
  if (colClustMethod == "pear_cor") {
    hc <- hclust(as.dist(1-cor(rmamat, method="pearson")), method="complete");
  } else if (colClustMethod == "euclidean") {
    hc <- hclust(dist(t(rmamat), method="euclidean"), method="complete");
  }
}

myheatcol <- colorpanel(40, "blue", "white", "red")

if (rowClust && colClust) {
  dendrogramArg = "both";
  hr <- as.dendrogram(hr);
  hc <- as.dendrogram(hc);
  permRmamat <- rmamat[,order.dendrogram(hc)][order.dendrogram(hr),];
} else if (rowClust && !colClust) {
  dendrogramArg = "row";
  hr <- as.dendrogram(hr)
  hc = FALSE
  permRmamat <- rmamat[order.dendrogram(hr),];
} else if (!rowClust && colClust) {
  dendrogramArg = "column"
  hr = FALSE
  hc <- as.dendrogram(hc)
  permRmamat <- rmamat[,order.dendrogram(hc)];
} else {
  dendrogramArg = "none"
  hr = FALSE
  hc = FALSE
  permRmamat <- rmamat;
}

if (numRows != 1){
  write.csv(permRmamat,outTxtFile)
} else {
  fixedMatrix = t(as.matrix(permRmamat))
  colnames(fixedMatrix) = names(permRmamat)
  rownames(fixedMatrix) = rownames(rma)
  write.csv(fixedMatrix,file=outTxtFile)
}

pdf(outFile, height = 2.5 + (0.5 * numRows), width = 3.6 + (0.1 * numCols), paper = "special", onefile = FALSE)

#if(bcolumn=='column'){
#    if (rowClustMethod != "none") {
#      heatmap.2(rmamat, Rowv=hr, Colv=hc, col=myheatcol, scale="column", density.info="none", trace="none", RowSideColors=mycolhc,dendrogram = dendrogramArg,cexCol=0.6,cexRow=0.6,margins=c(10,8))
#    } else {
#      heatmap.2(rmamat, Rowv=hr, Colv=hc, col=myheatcol, scale="column", density.info="none", trace="none",dendrogram = dendrogramArg,cexCol=0.6,cexRow=0.6,margins=c(10,8))
#    }
#} else{


if (bcolumn == "matrix") {
  legendText <- "Matrix Z-Score"
} else {
  legendText <- "Value"
}
if (numRows != 1 && numCols != 1){
  if (rowClustMethod != "none") {
        heatmap.2(rmamat, Rowv=hr,Colv=hc, col=myheatcol, scale="none", density.info="none", trace="none", RowSideColors=mycolhc,dendrogram = dendrogramArg,cexCol=0.6,cexRow=0.6,margins=c(10,8),key.xlab=legendText)
     
  } else {
        heatmap.2(rmamat, Rowv=hr,Colv=hc, col=myheatcol, scale="none", density.info="none", trace="none", dendrogram = dendrogramArg,cexCol=0.6,cexRow=0.6,margins=c(10,8),key.xlab=legendText)
  }
} else {
  image(rmamat)
}
#}
dev.off();
