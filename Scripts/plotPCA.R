#!/usr/bin/Rscript
# Usage: Rscript -i infile.covar -c component1-component2 -a annotation.file -o outfile.eps

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.5")

library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from ngsCovar)'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Components to plot'),
                    make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Annotation file is in plink cluster format

#################################################################################

# Read input file
covar <- read.table(opt$in_file, stringsAsFact=F);

# Read annot file
annot <- read.table(opt$annot_file, sep=" ", header=T); # note that plink cluster files are usually tab-separated instead

# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])



# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");


# Write eigenvalues
#write.table(eig, file = "eigen_scores.txt", quote = FALSE)


# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Tra <- factor(annot$IID)
PC$Lab <- factor(annot$FID)


# Write PC componen
write.table(PC, file = "PC_scores.txt", quote = FALSE)


title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop", shape="Tra")) + geom_text(data=PC, aes_string(x=x_axis, y=y_axis, label="Lab", vjust= -0.784), size = 2)+ ggtitle(title)
ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title)
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop", shape="Tra")) + scale_shape_manual(values = c(3,4,7,8,9,12,15,16,17)) + scale_color_manual(values=c("darkgreen", "red")) + ggtitle(title)
ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)

