### Rscript --vanilla pipeLine_DESeq2_comparison.R -c Control -t Treatment

## Notice
## Example visualization: The code below demonstrates the general structure of an RNA-versus-RPF correlation plot.
## The example code is simplified and does not represent the complete analysis implementation used in the study.


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
library(EnhancedVolcano)

setwd("./")

library("optparse")
 
option_list = list(
    make_option(c("-c", "--control"), type="character", default=NULL,
          help="control condition", metavar="character"),

    make_option(c("-t", "--treatment"), type="character", default=NULL,
          help="treatment condition", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

if (is.null(opt$control)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}



## pick up the selected samples
ctrl <- opt$control
treat <- opt$treatment

folder <- paste(treat,"_vs_",ctrl,sep="")
dir.create(file.path(paste(treat,"_vs_",ctrl,sep="")))

########  DE for RNA
count_table_RNA <- read.table("summary_gene_count_table.txt",head=T,check.names=F)
count_table_sample_all <- read.table("count_tab_sample_RNA.txt",head =T,sep="\t",check.names=F)
############# Step1 analysis RNA seq
count_table <- count_table_RNA
rownames(count_table) <- count_table$GeneID
count_table <- count_table[-1]
count_table_sample <- count_table_sample_all
row.names(count_table_sample)<-count_table_sample$sample_name
colnames(count_table)<-count_table_sample$sample_name

ctrl <- as.name(opt$control)
treat <- as.name(opt$treatment)

count_table_sample_2c <- count_table_sample[count_table_sample$cellLine %in% c(ctrl,treat),]
count_table_2c <- count_table[,rownames(count_table_sample_2c)]

dds <- DESeqDataSetFromMatrix(countData=count_table_2c, colData=count_table_sample_2c, design = ~ cellLine)
dim(dds)
dds <- dds[ rowMeans(counts(dds)) > 10, ]
dim(dds)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

res <- results(dds,contrast=c("cellLine",toString(treat),toString(ctrl)),independentFilter=FALSE)

resOrdered <- res[order(res$log2FoldChange),]
write.table(as.data.frame(resOrdered),file= paste(folder,"/diff_",treat,"_vs_",ctrl,"_RNA.txt",sep=""),row.names = T, sep="\t",quote=F)
diff_RNA <- as.data.frame(resOrdered)


#####################################################################################################################
#####################################################################################################################
####### RPF DE

count_table_RPF <- read.table("RPFs_table_geneLevel.txt",head=T,check.names=F)

count_table_sample_all <- read.table("count_tab_sample_RPF.txt",head =T,sep="\t",check.names=F)

############# Step1 analysis RPF seq
count_table <- count_table_RPF
count_table_sample <- count_table_sample_all
row.names(count_table_sample)<-count_table_sample$sample_name
colnames(count_table)<-count_table_sample$sample_name

ctrl <- as.name(opt$control)
treat <- as.name(opt$treatment)

count_table_sample_2c <- count_table_sample[count_table_sample$cellLine %in% c(ctrl,treat),]
count_table_2c <- count_table[,rownames(count_table_sample_2c)]

dds <- DESeqDataSetFromMatrix(countData=count_table_2c, colData=count_table_sample_2c, design = ~ cellLine)
dim(dds)
dds <- dds[ rowMeans(counts(dds)) > 10, ]
dim(dds)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

res <- results(dds,contrast=c("cellLine",toString(treat),toString(ctrl)),independentFilter=FALSE)
resOrdered <- res[order(res$log2FoldChange),]
write.table(as.data.frame(resOrdered),file= paste(folder, "/diff_",treat,"_vs_",ctrl,"_RPF.txt",sep=""),row.names = T, sep="\t",quote=F)
diff_RPF <- as.data.frame(resOrdered)


#######################################################################################################################
#################  functions
## Example visualization: The code below demonstrates the general structure of an RNA-versus-RPF correlation plot.
##The example code is simplified and does not represent the complete analysis implementation used in the study.

plot_RNA_RPF_correlation <- function(count_table, treat, ctrl) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install ggplot2")
  
  library(ggplot2)
  
  required <- c("log2FC_RNA", "log2FC_RPF")
  if (!all(required %in% names(count_table)))
    stop("Required columns are missing.")
  
  p <- ggplot(
    count_table,
    aes(x = log2FC_RNA, y = log2FC_RPF)
  ) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      x = paste0("RNA ", treat, " vs ", ctrl, " log2FC"),
      y = paste0("RPF ", treat, " vs ", ctrl, " log2FC"),
      title = "RNA vs RPF correlation"
    )
  
  return(p)
}

ID <- intersect(rownames(diff_RNA),rownames(diff_RPF))

table<- merge(as.data.frame(diff_RNA[ID,c(2,6)]), as.data.frame(diff_RPF[ID,c(2,6)]), by = "row.names", all = TRUE)
rownames(table) <- table$Row.names
table[is.na(table)] <- 0
newT <- table[ID,c(2:5)]

newName <- c("log2FC_RNA", "RNA_diff_padj", "log2FC_RPF", "RPF_diff_padj")
colnames(newT) <- newName

CorrO <- plot_RNA_RPF_correlation(newT,folder,treat,ctrl)
pdf(paste(folder,"/","DE_RPF vs DE_RNA"," (", treat,"_vs_", ctrl, ")_padj.pdf", sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
CorrO
garbage <- dev.off()


########################################################################################*****************************
########################################################################################*****************************
## based on p-value
##########################################################################################################################

ID <- intersect(rownames(diff_RNA),rownames(diff_RPF))

table<- merge(as.data.frame(diff_RNA[ID,c(2,5)]), as.data.frame(diff_RPF[ID,c(2,5)]), by = "row.names", all = TRUE)
rownames(table) <- table$Row.names
table[is.na(table)] <- 0
newT <- table[ID,c(2:5)]
newName <- c("log2FC_RNA", "RNA_diff_pvalue", "log2FC_RPF", "RPF_diff_pvalue")
colnames(newT) <- newName

CorrO <- plot_RNA_RPF_correlation (newT,folder,treat,ctrl)
pdf(paste(folder,"/","DE_RPF vs DE_RNA"," (", treat,"_vs_", ctrl, ")_pvalue.pdf", sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
CorrO
garbage <- dev.off()
