
## Exploratory Data Analysis

rm(list=ls())

#Do an exploratory data analysis of a matrix of expression values.

## Data Import

anno = read.table("SampleAnnotation.txt", as.is=TRUE, sep="\t", quote="", row.names=1, header=TRUE)

x = read.table("expressiondata.txt", as.is=TRUE, sep="\t", quote="", row.names=1, header=TRUE, check.names = FALSE)
x = as.matrix(x)


## Define samples and colors and phenotype

source("https://bioconductor.org/biocLite.R")
biocLite("limma")

str(anno)
ysamples = rownames(anno)  #gives names of rows
#ysamples
#nrow(anno)  #gives number of rows
colors = rainbow(nrow(anno))

isNorm = anno$TissueType == "norm"
#isNorm
isSick = anno$TissueType == "sick"
isAcute = anno$TissueType == "acute"


#  * distributions: *boxplot*, *density*, *limma::plotDensities*

par(mfrow=c(3,3))
boxplot(log2(x), main="Expressionvalues: Boxplot", col=colors, ylab = "log2 of expressionvalues",
        xlab = "Sample names", las=2)  #take log 2 of matrix before boxplot, las=2 sets axis-labels to perpendicular

#?boxplot
#boxplot(x[1,]~x,data=x, main="Boxplot",
#   xlab="Samples", ylab="Expression values")

plot(density(log2(x)),  main="Expressionvalues: Density", ylab = "Density of expressionvalues",
     xlab = "log2 of expressionvalues", las=2)

limma::plotDensities(log2(x), legend = "topright")
#legend(location="topleft")


#  * normalization: *limma::normalizeQuantiles*

xnorm <- limma::normalizeQuantiles(x)

boxplot(log2(xnorm), main="Normalized Expressionvalues: Boxplot", col=colors, ylab = "Epressionvalues",
        xlab = "Sample names", las=2)

plot(density(log2(xnorm)),  main="Normalized Expressionvalues: Density", ylab = "Density of expressionvalues",
     xlab = "Expressionvalues", las=2)

limma::plotDensities(log2(xnorm), legend = "topright")


#  * clustering: *hclust*

#?hclust
#?dist
#head(x)
#x[1,]
clust <- hclust(dist(x[1,]),"cen")  #clustering of first gene
#?cor
#?as.dist
#clust <- hclust(dist((1-cor(x))/2))
#plot(clust)
clust <- hclust(as.dist((1-cor(x))/2))
plot(clust)
#all the norms are grouped (with one exception of a sick) and all the acutes and sicks are clustered.

#  * heatmap: *heatmap.2* or *pheatmap*

#?heatmap.2
heatmap.2(log2(x))
#?pheatmap
pheatmap(log2(x))

#  * correlation matrix: *cor* and *image*

#?cor
#?image
correl <- cor(x)
image(correl, main="visual representation of the correlation matrix", col = colors)

#  * reduced dimensionality representation: *cmdscale* and *prcomp*

#?cmdscale
cmdsc <- cmdscale(as.dist((1-cor(x))/2), k=17)
plot(cmdsc)

#?prcomp
pca <- prcomp(log2(x))
plot(pca)
