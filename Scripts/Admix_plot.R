
args <- commandArgs(TRUE)
infile <- args[1]
info <- args[2]
n <- args[3]

pop<-read.table(info,as.is=T)
admix<-t(as.matrix(read.table(infile)))
admix<-admix[,order(pop[,1])]
pop<-pop[order(pop[,1]),]
h<-barplot(admix,col=rainbow(n),space=0,border=NA,xlab="Individuals",ylab="admixture")
text(tapply(1:nrow(pop),pop[,1],mean), -0.07, unique(pop[,1]),xpd=T, srt = -90)
