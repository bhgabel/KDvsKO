library(tidyverse)
library("DESeq2")

#Raw RNA-seq data
kd <- read.csv("KD_RawReadCounts_RNAseq_MCF7.csv", header=T)
gene_info <- select(kd, "geneID", "GeneName")

kd <- select(kd, "geneID", "shLuc.Veh.1", "shLuc.veh.2",
             "shMLH1.veh.1", "shMLH1.veh.2")
kd <- column_to_rownames(kd, var="geneID")


#metadata for DESeq
metadata <- data.frame(condition=c("shLuc", "shLuc", "shMLH1", "shMLH1"))
rownames(metadata) <- colnames(kd)

all(colnames(kd)[-1] %in% rownames(metadata))

#create DESeq object
dds <- DESeqDataSetFromMatrix(countData=kd, colData=metadata, design=~condition)

# filter any counts less than 10
dds <- dds[rowSums(counts(dds)) >= 10,]

#Differential expression
#defaults to alpha=0.01
dds <- DESeq(dds)
res <- results(dds)

summary(results)

#count of genes with p-value < 0.01
sum(res$padj < 0.1, na.rm=TRUE)

plotPCA(rlog(dds))

results <- results %>%
  mutate(Expression = case_when(log2FoldChange >= log(1) & padj <= 0.01 ~ "Up-regulated",
                           log2FoldChange <= -log(1) & padj <= 0.01 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))

output <- rownames_to_column(results, var="geneID")
output <- left_join(output, gene_info, by="geneID")
write.csv(subset(output, padj <= 0.01), "KD_DESeq_results.csv", row.names=F)

#Volcano plot

top <- 10
# we are getting the top 10 up and down regulated genes by filtering the column Up-regulated and Down-regulated and sorting by the adjusted p-value. 
top_genes <- bind_rows(
  results %>%
    filter(Expression == 'Up-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top),
  results %>%
    filter(Expression == 'Down-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top)
)
# create a datframe just holding the top 10 genes
Top_Hits <-  head(arrange(results,pvalue),10)
Top_Hits


# basic plot with line + red for p < 0.05
ggplot(results, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  #geom_hline(yintercept= -log10(0.01), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  theme(legend.position = "none")
  geom_text_repel(aes(label = label), size = 2.5)
