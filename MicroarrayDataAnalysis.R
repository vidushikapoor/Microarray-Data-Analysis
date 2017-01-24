#Expression of survivin in lung eosinophils is associated with pathology in a mouse model of allergic asthma
#Model: BALB

#Microarray Data Normalization
#Bioconductor Package
source("http://bioconductor.org/biocLite.R") 
biocLite("marray")
library(marray)

#Set working directory
setwd("/Users/vidushikapoor/Documents/3 sem Bioinfo/Statistics/project")

#Read the target file 
MouseTargets <- read.marrayInfo("MouseTargetsBAOP.txt",info.id = "FileName")
mouse.names <- read.marrayInfo("GPL3069-tbl-1.txt", info.id = 6:7, label = 7, skip = 0)

#Read the .spot files
mraw <- read.Spot(targets = MouseTargets)

#Produce quality assessment plots
jpeg('boxplot-normc57-control.jpg')
boxplot(mraw)
dev.off()

#Background images to see artifacts
image(mraw[,1],xvar="maRb")

#Conduct Normalization

#Global Median Normalization
mraw.m.norm<-maNorm(mraw,norm="m")
jpeg('boxplot+Median+normc57+control.jpg')
boxplot(mraw.m.norm)
dev.off()

#Global Loess normalization
mraw.l.norm <- maNorm(mraw, norm="l")
jpeg('boxplot+Loess+normc57+control.jpg')
boxplot(mraw.l.norm)
dev.off()

#Print-Tip Loess Normalization
mraw.p.norm <- maNorm(mraw, norm="p")
jpeg('boxplot+Tip+normc57+control.jpg')
boxplot(mraw.p.norm)
dev.off()

#Look up the normalized values
maM(mraw.m.norm)
head(maM(mraw.m.norm))

# Limma package

source("http://bioconductor.org/biocLite.R")
biocLite("limma",dependencies=TRUE)
biocLite("convert")
library(limma)
library(convert)

spot.limma<-as(mraw,"RGList")

#Normalization in Limma
MA<-normalizeWithinArrays(spot.limma)
class(MA)
MA<-normalizeBetweenArrays(MA)

#Specify design matrix
design<-c(1,1,-1,-1)
design

#Fitting linear model
fit<-lmFit(MA,design)
class(fit)
ordinary.t<-fit$coef/fit$stdev.unscaled / fit$sigma

#Shrink gene-wise standard deviations and conduct t-tests
fit<-eBayes(fit)

#Top table
a<-topTable(fit, number = length(fit), adjust.method = "BH")
foo<-list()
temp=0.05
for(i in 1:length(fit)){
  if (a$P.Value[i]<temp)
  {
  foo<-append(foo,a$Name[i])
  }
}

lapply(foo, write, file="PBS22.text", append=T)
write(a$Name, file = "data.txt")

#Volcano plot
jpeg('volcanoPlotc57-control.jpg')
volcanoplot(fit, highlight = 30, names = fit$genes$Name)
dev.off()

#Annotated MA plot
ord <- order(fit$lods, decreasing = TRUE)
top30 <- ord[1:30]
top3text(fit$Amean[top30], fit$coef[top30], labels = fit$genes[top30, "Name"], cex = 0.8, col = "blue")

