#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

F1 = args[1]
F2 = paste(F1, ".pdf", sep="")

#function to normalize
norm <- function(x) x/sum(x)
#read data
sfs <- (scan(F1))
#the variability as percentile
pvar<- (1-sfs[1]-sfs[length(sfs)])*100
#the variable categories of the sfs
sfs<-norm(sfs[-c(1,length(sfs))]) 
pdf(F2)
barplot(sfs,legend=paste("Variability:= ",round(pvar,3),"%"),xlab="number of derived sites",
names=1:length(sfs),ylab="Proportions",main=F1,col='brown')



