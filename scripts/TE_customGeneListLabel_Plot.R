
### Rscript --vanilla TE_customGeneListPlot.R -c Control -t Treatment
#### Note: need to be modified using variable label same as input to specify custom condition names in above command line


library("RColorBrewer")
library("ggplot2")
library("plyr")
library("NMF")
library(ggplot2)
library(tibble)
library(magrittr)
library(stringr)
library("ggrepel")
library(DESeq2)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)

setwd("./")

library("optparse")

# example only
SpecificGeneList <- c(
  "Myc",
  "Cdk7"
)


option_list = list(
    make_option(c("-c", "--control"), type="character", default=NULL,
          help="control condition", metavar="character"),

    make_option(c("-t", "--Treatment"), type="character", default=NULL,
          help="Treatment condition", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

if (is.null(opt$control)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


ctrl <- opt$control

treat <- opt$Treatment

#### Note: the following condition labels needs to be modified to variable label as above for pipeline ready run on the specific project
####

dir.create(file.path(paste(treat,"_vs_",ctrl,sep="")))
folder <- paste(treat,"_vs_",ctrl,sep="")

obj <- readRDS("TE_analysis_full.rds")  ## please refer to the TE analysis result obj

dds <- obj$dds
res <- obj$res
resLFC <- obj$resLFC

# remove NA values
res <- res[!is.na(res$stat), ]


################### volcano plot
pdf(paste(folder, "/TE_Treatment_vs_Ctrl_VolcanoPlot.pdf", sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
# Define colors
col_vec <- ifelse(
  res$log2FoldChange > 0 & res$padj < 0.05, "red",        # upregulated significant
  ifelse(res$log2FoldChange < 0 & res$padj < 0.05, "blue",  # downregulated significant
         "grey")                                            # not significant
)
# Replace NA padj safely
padj_safe <- res$padj
padj_safe[is.na(padj_safe)] <- 1
plot(
  res$log2FoldChange,
  -log10(padj_safe),
  pch=20,
  cex=0.6,
  col=col_vec,
  xlab="log2 TEs fold change (Treatment_TE vs Ctrl_TE)",
  ylab="-log10 adjusted P-value",
  main="Treatment_TE vs Ctrl_TE Volcano Plot"
)
abline(v=0, col="black", lty=2)
abline(h=-log10(0.05), col="blue", lty=2)
dev.off()


############## volcano plot with specific genes labeling
# contain for the specific study


