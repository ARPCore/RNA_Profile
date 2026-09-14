### Note: need to introduce the corresponding data/results from DE/TE analysis

library("RColorBrewer")
library("ggplot2")
library("plyr")
library("NMF")
library(ggplot2)
library(tibble)
library(magrittr)
library(stringr)
library("ggrepel") 

library("ReactomePA")
library("DOSE")
library(fgsea)
library(clusterProfiler)
library(clusterProfiler)
library(org.Mm.eg.db)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)

geneList <- res$stat
names(geneList) <- rownames(res)

# sort decreasing
geneList <- sort(geneList, decreasing = TRUE)

gene_df <- bitr(names(geneList),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

#Then rebuild geneList:
geneList <- geneList[gene_df$SYMBOL]
names(geneList) <- gene_df$ENTREZID

#Step 3 — Remove duplicates & sort
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

cgsea_res <- gseGO(geneList = geneList,
                     ont = "BP",
                     OrgDb = org.Mm.eg.db,
                     minGSSize = 15,
                     maxGSSize = 500,
                     eps = 0,
                     nPermSimple = 10000,
                     seed = TRUE)

write.table(cgsea_res@result,file="GO_results/KOvsWT_GSEA_GOterms.txt",sep="\t",quote=F) ### except immune in wpTerm list

#########################################################################################################
#########################################################################################################
df <- cgsea_res@result %>% arrange(p.adjust) %>% head(12)
pdf("GO_results/GSEA_GO_DotPlot.pdf")
par(mar=c(5.1,4.1,4.1,4.1))
ggplot(df, aes(x = NES,
               y = reorder(Description, NES))) +
  geom_point(aes(size = setSize,
                 color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "Top 12 Enriched GO Biological Processes",
    x = "Normalized Enrichment Score (NES)",
    y = "",
    color = "Adjusted p-value",
    size = "Gene Set Size"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()


##########################################################################################
### gene set enrichment plot on
library(enrichplot)
library(ggplot2)

p <- gseaplot2(
cgsea_res,
  geneSetID = 1:5,
  title = NULL
)

pdf("GO_results/Gene_set_enrichment_plot_multiple.pdf")
par(mar=c(5.1,4.1,4.1,4.1))
print(p)
dev.off()
