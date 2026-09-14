
### Rscript --vanilla TE_customGeneListPlot.R -c Control -t Treatment

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

df <- data.frame(
  gene = rownames(res),
  log2FC = res$log2FoldChange,
  padj = res$padj
)

df$negLog10 <- -log10(df$padj)

df$label <- ifelse(df$gene %in% SpecificGeneList, df$gene, NA)

pdf(paste(folder, "/TE_Treatment_vs_Ctrl_VolcanoPlot_SpecificGenesLabeling.pdf", sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
padj_safe <- res$padj
padj_safe[is.na(padj_safe)] <- 1
yvals <- -log10(padj_safe)
yvals[is.infinite(yvals)] <- max(yvals[is.finite(yvals)], na.rm=TRUE)
col_vec <- ifelse(
  res$log2FoldChange > 0 & res$padj < 0.05, "red",
  ifelse(res$log2FoldChange < 0 & res$padj < 0.05, "blue", "grey")
)
plot(
  res$log2FoldChange,
  yvals,
  pch=20,
  cex=0.6,
  col=col_vec,
  xlab="log2 TEs fold change (Treatment_TE vs Ctrl_TE)",
  ylab="-log10 adjusted P-value",
  main="Treatment_TE vs Ctrl_TE Volcano Plot with Gene Labels"
)
abline(v=0, col="black", lty=2)
abline(h=-log10(0.05), col="blue", lty=2)
# label selected genes
idx <- which(rownames(res) %in% SpecificGeneList)
text(
  x = res$log2FoldChange[idx],
  y = yvals[idx],
  labels = rownames(res)[idx],
  pos = 3,
  cex = 0.8,
  col = "black"
)
dev.off()


#################################### remove NA padj for plot
df <- data.frame(
  gene = rownames(res),
  log2FC = res$log2FoldChange,
  padj = res$padj
)
# flag NA before filtering (important!)
df$padj_status <- ifelse(is.na(df$padj), "NA", "OK")
# remove NA for volcano stats
df_plot <- df[!is.na(df$padj), ]
df_plot$negLog10 <- -log10(df_plot$padj)
df_plot$label <- ifelse(df_plot$gene %in% SpecificGeneList, df_plot$gene, NA)

df_plot$color <- ifelse(
  df_plot$log2FC > 0 & df_plot$padj < 0.05, "Up in Treatment_TE",
  ifelse(df_plot$log2FC < 0 & df_plot$padj < 0.05, "Down in Treatment_TE", "Not Significant")
)

pdf(paste(folder, "/TE_Treatment_vs_Ctrl_VolcanoPlot_SpecificGenesLabeling_ggrepel.pdf", sep=""))
ggplot(df_plot, aes(x = log2FC, y = negLog10)) +
  geom_point(aes(color = color), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up in Treatment_TE"="red", "Down in Treatment_TE"="blue", "Not Significant"="grey")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  theme_classic() +
  labs(
    x = "log2 TE fold change (Treatment_TE vs Ctrl_TE)",
    y = "-log10 adjusted P-value",
    title = "TEs Volcano Plot (Treatment_TE vs Ctrl_TE)"
  )

dev.off()


