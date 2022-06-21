
fin <- commandArgs(T)
fout <- paste(fin,".pdf",sep="",collaps="")

library(ape)
library(phangorn)

trees <- read.tree(fin, skip=2)
tree <- read.tree(fin)[[1]]

pdf(file=fout)

cons = plotBS(tree, trees, type="phylo", cex=0.5)

invisible(dev.off())

write.tree(cons, file="noccaea_21_boot_cons.tre")

