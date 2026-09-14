### Rscript --vanilla pipeline_TE_byConditions.R -c Control -t Treatment
# note: certain features which were specifically used in the study were removed


library("RColorBrewer")
library("ggplot2")
library("plyr")
library("NMF")
library(ggplot2)
library(tibble)
library(magrittr)
library(stringr)
library("ggrepel") ######### Provides text and label geoms for 'ggplot2' that help to avoid
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)

setwd("./")

library("optparse")
 
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


Ctrl <- opt$control

Treatment <- opt$Treatment

dir.create(file.path(paste(Treatment,"_vs_",Ctrl,sep="")))
folder <- paste(Treatment,"_vs_",Ctrl,sep="")

## Read RNA table
count_table_RNA <- read.table("summary_gene_count_table.txt",head=T,check.names=F)
count_table_sample_all <- read.table("count_tab_sample_RNA.txt",head =T,sep="\t",check.names=F)
count_table <- count_table_RNA
rownames(count_table) <- count_table$GeneID
count_table <- count_table[-1]
count_table_sample <- count_table_sample_all
row.names(count_table_sample)<-count_table_sample$sample_name
colnames(count_table)<-count_table_sample$sample_name
RNA_meta <- count_table_sample
RNA_count <- count_table

## Read RPF table
count_table_RNA <- read.table("RPFs_table_geneLevel.txt",head=T,check.names=F)
count_table_sample_all <- read.table("count_tab_sample_RPF.txt",head =T,sep="\t",check.names=F)
count_table <- count_table_RNA
count_table_sample <- count_table_sample_all
row.names(count_table_sample)<-count_table_sample$sample_name
colnames(count_table)<-count_table_sample$sample_name
RPF_meta <- count_table_sample
RPF_count <- count_table

## REMOVE sum of 0 for RPF and RNA
RPF_count <- RPF_count[rowSums(RPF_count)>0,]
RNA_count <- RNA_count[rowSums(RNA_count)>0,]

table <- merge(RPF_count,RNA_count, by ="row.names", all = TRUE)
table <- table[complete.cases(table), ] ## remove NA
rownames(table) <- table$Row.names

############################################
## function to get splited sample condition
spl_name <- function(x) unlist(str_split(x,"_"))[2]
table_names <- names(table)
table_Ctrl <- table[,grepl(Ctrl, table_names)]
table_Ctrl_RPF <- table_Ctrl[,grepl("RPF", names(table_Ctrl))]
table_Ctrl_RNA <- table_Ctrl[,grepl("RNA", names(table_Ctrl))]

# "Treatment" sample to compare to
table_Treatment <- table[,grepl(Treatment, table_names)]
table_Treatment_RPF <- table_Treatment[,grepl("RPF", names(table_Treatment))]
table_Treatment_RNA <- table_Treatment[,grepl("RNA", names(table_Treatment))]


##########################################################################################################################
########### Calculate TE for Ctrl RPF vs RNA
Ctrl_folder <- paste(folder,"/TE_Condition_",Ctrl,sep="")
dir.create(file.path(Ctrl_folder))

Ctrl_S <- cbind(table_Ctrl_RPF,table_Ctrl_RNA)
Ctrl_S <- Ctrl_S[rowSums(Ctrl_S)>0,]

tt <- rbind(RPF_meta,RNA_meta)
meta_names <- rownames(tt)
meta_Ctrl <- tt[grepl(Ctrl,meta_names),]

########## all data
dds <- DESeqDataSetFromMatrix(countData=Ctrl_S, colData=meta_Ctrl, design = ~seqType)
dds <- dds[ rowMeans(counts(dds)) > 10, ]
dim(dds)
dds <- DESeq(dds)
nor_counts <- counts(dds,normalized=T)

tt <- as.data.frame(nor_counts)
dim(tt)

############################################
write.table(tt ,file=paste(Ctrl_folder,"/",Ctrl,"_RPFvsRNA_geneLevel_ddsNormalization.txt",sep=""),row.names=T,sep="\t", quote=F)

######## Generate Volcano plot
res <- results(dds,contrast=c("seqType","RPF","RNA"),independentFilter=FALSE)
resOrdered <- res[order(res$log2FoldChange),]

write.table(resOrdered ,file=paste(Ctrl_folder,"/diffGene_",Ctrl,"_RPFvsRNA_all.txt",sep=""),row.names = T, sep="\t",quote=F)

### out put diff table for Ctrl sample
diff_table_Ctrl <- resOrdered

################## customized
## create custom key-value pairs for 'pj>=.05', 'pf<.05$|log2FC|<=2', 'pf<.05$log2FC>2' pf<.05$log2FC<-2
# this can be achieved with nested ifelse statements
res <- na.omit(res) ## remove NA in pvalue or adjpvalue

C1 <- res[res$padj>=0.05,]
C2 <- res[(res$padj<0.05 & abs(res$log2FoldChange)<=log2(1.5)),]
C3 <- res[(res$padj<0.05 & res$log2FoldChange>log2(1.5)),]
C4 <- res[(res$padj<0.05 & res$log2FoldChange< -log2(1.5)),]

write.table(C3,file=paste(Ctrl_folder,"/diffGene_",Ctrl,"_RPFvsRNA_UpregulateGenes_ByAdjP.txt",sep=""),row.names = T, sep="\t",quote=F)
write.table(C4,file=paste(Ctrl_folder,"/diffGene_",Ctrl,"_RPFvsRNA_DownregulateGenes_ByAdjP.txt",sep=""),row.names = T, sep="\t",quote=F)

####### include all
pdf(paste(Ctrl_folder,"/TE_VolcanoPlot_",Ctrl,"_RPFvsRNA_GeneLevel_basedOnAdjustedPValue_FC_1ndHalf_labeled.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    #selectLab = c("Hist2h3c2","Crabp1"),
    xlab = bquote(~Log[2]~ "TE"),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Ctrl,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    #xlim = c(-15,15),
    labSize = 2.0,#3
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10, or 6
    legendIconSize = 3, ## 3
    #DrawConnectors = FALSE,
    #widthConnectors = 0.5,
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
)
garbage <- dev.off()

pdf(paste(Ctrl_folder,"/TE_VolcanoPlot_",Ctrl,"_RPFvsRNA_GeneLevel_basedOnAdjustedPValue_FC_1ndHalf.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = NA,
    x = "log2FoldChange",
    y = "padj",
    #selectLab = c("Hist2h3c2","Crabp1"),
    xlab = bquote(~Log[2]~ "TE"),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Ctrl,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    #xlim = c(-15,15),
    labSize = 2.0,#3
    colAlpha = 1,
    #legend=c("NS","Log2 FC","Adjusted p-value",
    #        "Adjusted p-value & Log2 FC"),
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10, or 6
    legendIconSize = 3, ## 3
    #DrawConnectors = FALSE,
    #widthConnectors = 0.5,
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    
    #colCustom = keyvals
)
garbage <- dev.off()


##################
################## Based on p-Value
C1 <- res[res$pvalue>=0.05,]
C2 <- res[(res$pvalue<0.05 & abs(res$log2FoldChange)<=log2(1.5)),]
C3 <- res[(res$pvalue<0.05 & res$log2FoldChange>log2(1.5)),]
C4 <- res[(res$pvalue<0.05 & res$log2FoldChange< -log2(1.5)),]

write.table(C3,file=paste(Ctrl_folder,"/diffgene_",Ctrl,"_RPFvsRNA_UpregulateGenes_byPvalue.txt",sep=""),row.names = T, sep="\t",quote=F)
write.table(C4,file=paste(Ctrl_folder,"/diffgene_",Ctrl,"_RPFvsRNA_DownregulateGenes_byPvalue.txt",sep=""),row.names = T, sep="\t",quote=F)

####### include all
pdf(paste(Ctrl_folder,"/TE_VolcanoPlot_",Ctrl,"_RPFvsRNA_GeneLevel_basedOnPValue_FC_1ndHalf_labeled.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~italic(P-value)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Ctrl,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    labSize = 2.0,#3
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10
    legendIconSize = 3, ## 3
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    
    #colCustom = keyvals
)
garbage <- dev.off()

pdf(paste(Ctrl_folder,"/TE_VolcanoPlot_",Ctrl,"_RPFvsRNA_GeneLevel_basedOnPValue_FC_1ndHalf.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = NA,
    x = "log2FoldChange",
    y = "pvalue",
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~italic(P-value)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Ctrl,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    #xlim = c(-15,15),
    labSize = 2.0,#3
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10
    legendIconSize = 3, ## 3
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
)
garbage <- dev.off()


########################################################################################*****************************
########################################################################################*****************************
## "Treatment" sample to compare to
##########################################################################################################################
########### Calculate TE for Treatment RPF vs RNA
# Prepare tables
## create separate table
Treatment_folder <- paste(folder,"/TE_Condition_",Treatment,sep="")
dir.create(file.path(Treatment_folder))

Treatment_S <- cbind(table_Treatment_RPF,table_Treatment_RNA)
Treatment_S <- Treatment_S[rowSums(Treatment_S)>0,]

tt <- rbind(RPF_meta,RNA_meta)
meta_names <- rownames(tt)
meta_Treatment <- tt[grepl(Treatment,meta_names),]

########## all data
dds <- DESeqDataSetFromMatrix(countData=Treatment_S, colData=meta_Treatment, design = ~seqType)
dds <- dds[ rowMeans(counts(dds)) > 10, ]
dim(dds)
dds <- DESeq(dds)
nor_counts <- counts(dds,normalized=T)

tt <- as.data.frame(nor_counts)
dim(tt)

############################################
write.table(tt ,file=paste(Treatment_folder,"/",Treatment,"_RPFvsRNA_geneLevel_ddsNormalization.txt",sep=""),row.names=T,sep="\t", quote=F)

######## Generate Volcano plot
res <- results(dds,contrast=c("seqType","RPF","RNA"),independentFilter=FALSE)
resOrdered <- res[order(res$log2FoldChange),]

write.table(resOrdered ,file=paste(Treatment_folder,"/diffGene_",Treatment,"_RPFvsRNA_all.txt",sep=""),row.names = T, sep="\t",quote=F)

diff_table_Treatment <- resOrdered

################## customized
res <- na.omit(res) ## remove NA in pvalue or adjpvalue
####### include all

pdf(paste(Treatment_folder,"/TE_VolcanoPlot_",Treatment,"_RPFvsRNA_GeneLevel_basedOnAdjustedPValue_FC_1ndHalf.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = NA,
    x = "log2FoldChange",
    y = "padj",
    xlab = bquote(~Log[2]~ "TE"),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Treatment,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    labSize = 2.0,#3
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10, or 6
    legendIconSize = 3, ## 3
    #DrawConnectors = FALSE,
    #widthConnectors = 0.5,
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    
)
garbage <- dev.off()

##################
################## Based on p-Value

pdf(paste(Treatment_folder,"/TE_VolcanoPlot_",Treatment,"_RPFvsRNA_GeneLevel_basedOnPValue_FC_1ndHalf.pdf",sep=""))
par(mar=c(5.1,4.1,8.1,4.1))
EnhancedVolcano(res,
    lab = NA,
    x = "log2FoldChange",
    y = "pvalue",
    xlab = bquote(~Log[2]~ "TE"),
    ylab = bquote(~-Log[10]~italic(P-value)),
    axisLabSize = 12,
    title = paste("VolcanoPlot_",Treatment,"_gene_RPFvsRNA",sep=""),
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 1,#
    labSize = 2.0,#3
    colAlpha = 1,
    legendPosition = "bottom",
    legendLabSize = 4.5,## 10
    legendIconSize = 3, ## 3
    border = "full",
    borderWidth = 1,
    borderColour = "black",
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    
    colCustom = keyvals
)
garbage <- dev.off()


########################################################################################*****************************
########################################################################################*****************************
## TEs interactive analysis
##########################################################################################################################
###########

dir.create(file.path(paste(folder,"/","TEs_Correlation",sep="")))
folder_te <- paste(folder,"/","TEs_Correlation",sep="")

########################################################################################*****************************
########################################################################################*****************************
## TE correlation analysis based on the DEs based on Fold Change
correlationFigureByColor <- function(count_table) {
    aa<-count_table[,c(1,3)]
    x_name <- names(aa)[1]
    y_name <- names(aa)[2]
    aa$Treatment_TEvsCtrl_TE <- aa[,2]-aa[,1]
    aa$Category <- ifelse(aa$Treatment_TEvsCtrl_TE> 1, paste(Treatment,"_TEvs",Ctrl,"_TE>1",sep=""),
                    ifelse(aa$Treatment_TEvsCtrl_TE< -1, paste(Treatment,"_TEvs",Ctrl,"_TE< -1",sep=""), paste("|",Treatment,"_TEvs",Ctrl,"_TE|<=1",sep="")))
    ## Export TE with significance table
    write.table(aa,file=paste(folder_te,"/TE_Category_Compare_caculation_table.txt",sep=""),row.names=T,sep="\t")
    fit <- lm(aa[,2]~aa[,1]-1, data=aa)
    plot <- ggplot(aa, aes_string(x = colnames(aa)[1], y = colnames(aa)[2])) +
            geom_point(aes(color = Category)) +
            labs(color = paste(Treatment,"_TEvs",Ctrl,"_TE (log2)",sep=""),title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
            " Slope =",signif(fit$coef[[1]], 5)))+
            scale_color_manual(values = c("grey", "red","blue")) +
            theme_bw(base_size = 10) + theme(legend.position = "bottom") +
            xlab(expr(paste(Delta,!!Ctrl,"_log2(RPF_avg/RNA_avg)",sep=""))) + ylab(expr(paste(Delta,!!Treatment,"_log2(RPF_avg/RNA_avg)",sep="")))+
            geom_abline(intercept=0, slope=fit$coefficients[1], color="blue",size=0.6,linetype="twodash") +
            geom_abline(intercept = 0, slope = 1,color="black",size=0.6) +
            expand_limits(x = 0, y = 0)
            
    return(plot)
}


ID <- intersect(rownames(diff_table_ctrl),rownames(diff_table_treat))
table<- merge(as.data.frame(diff_table_ctrl[ID,c(2,6)]), as.data.frame(diff_table_treat[ID,c(2,6)]), by = "row.names", all = TRUE)
rownames(table) <- table$Row.names
table[is.na(table)] <- 0
newT <- table[ID,c(2:5)]
newName <- c(paste("log2FC_",ctrl,"_RPFvsRNA",sep=""), paste(ctrl,"_padj",sep=""),paste("log2FC_",treat,"_RPFvsRNA",sep=""), paste(treat,"_padj",sep=""))
colnames(newT) <- newName

write.table(newT ,file=paste(folder_te,"/", treat, "_vs_", ctrl, "_TE_Calculation_table_RPFvsRNA.txt",sep=""),row.names=T,sep="\t", quote=F)


CorrO_FC <- correlationFigureByColor(newT)
pdf(paste(folder_te,"/Delta_",Treatment," vs Delta_", Ctrl, " (", Treatment, "_TEvs",Ctrl,"_TE).pdf",sep=""))
par(mar=c(5.1,4.1,4.1,4.1))
CorrO_FC
garbage <- dev.off()

########################################################################################*****************************
########################################################################################*****************************
## TE correlatiopn analysis based on the TE of sample based conditions
##########################################################################################################################

T_table <- cbind(table_Ctrl, table_Treatment)
sample_info <- rbind(meta_Ctrl, meta_Treatment)

## Relevel before DESeq()
sample_info$condition <- relevel(factor(sample_info$condition), ref = "Control")
dds <- DESeqDataSetFromMatrix(countData = T_table, colData   = sample_info, design    = ~ condition + seqType + condition:seqType)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="conditionTreatment.seqTypeRPF")
resOrdered <- res[order(res$log2FoldChange),]
write.table(resOrdered, file=paste(folder_te,"/", Treatment, "_vs_", Ctrl, "_TEs_Interaction_table.txt",sep=""),row.names = T, sep="\t",quote=F)

# Optional: shrink LFCs for visualization
resLFC <- lfcShrink(dds, coef="conditionTreatment.seqTypeRPF", type="ashr")  # needs ashr package

### Important:
saveRDS(list(dds = dds,res = res,resOrdered = resOrdered, resLFC = resLFC), file = paste0(folder_te, "/TE_analysis_full.rds"))
## provide the TE analysis Object for downstream TEs comparison analysis
