library(DESeq2)
library(tximport)
library(yaml)
library(dplyr)
library(biomaRt)
library(ggfortify)
library(RColorBrewer)

# Read metadata
metadata <- read.table("../metadata/P2022_metadata.tsv", sep="\t", header = T)
config <- unlist(read_yaml("../config.yaml"))

kal_dirs <- dir(file.path(config["data.kallisto"]))
sample_id <- sapply(strsplit(kal_dirs,"_"), `[`, 1)
sample_id <- sapply(strsplit(sample_id,"-"), `[`, 1)
kal_dirs <- paste(config["data.kallisto"], kal_dirs, sep="/")
kal_df <- data.frame(KallistoDirs = kal_dirs, SampleID = sample_id)
names(kal_dirs) <- metadata$SampleID

metadata <- merge(metadata, kal_df, by="SampleID")

# Create transcript to gene mapping

ensembl <- useEnsembl(biomart = "genes", version=109)
ensembl <- useDataset(dataset = "rnorvegicus_gene_ensembl", mart = ensembl)
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "rgd_symbol", "entrezgene_description"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, symbol = rgd_symbol, description=entrezgene_description)

# Import kallisto data
files <- file.path(kal_dirs, "abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", tx2gene=t2g[,1:2], ignoreTxVersion = T, txOut = F)
head(txi.kallisto$counts)

# PCA plots
pca <- prcomp(t(log2(txi.kallisto$counts+0.5)))
autoplot(pca, data = metadata, colour = 'Group')
ggsave("../plots/PCA_all.png", width=7, height=4)

ind <- grep("Auge", metadata$Group, invert=T)
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_no_auge.png", width=7, height=6)

ind <- which(!grepl("Auge", metadata$Group) & grepl("np", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_np_no_auge.png", width=5, height=3)

ind <- which(!grepl("Auge", metadata$Group) & grepl("pp", metadata$Group) & grepl("SD", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_pp_SD_no_auge.png", width=5, height=3)

ind <- which(!grepl("Auge", metadata$Group) & grepl("d21", metadata$Group) & grepl("SD", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_d21_SD_no_auge.png", width=5, height=3)

ind <- which(!grepl("Auge", metadata$Group) & grepl("pp", metadata$Group) & grepl("PE", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_pp_PE_no_auge.png", width=5, height=3)
# Remove outlier
ind <- ind[-which(ind==105)]
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_pp_PE_no_auge_no_outlier.png", width=5, height=3)

ind <- which(!grepl("Auge", metadata$Group) & grepl("d21", metadata$Group) & grepl("PE", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_d21_PE_no_auge.png", width=5, height=3)

ind <- which(grepl("LV", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_LV.png", width=5, height=3)

ind <- which(grepl("RV", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_RV.png", width=5, height=3)

ind <- which(grepl("Ap", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_Ap.png", width=5, height=3)

ind <- which(grepl("Sept", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_Sept.png", width=5, height=3)

ind <- which(grepl("Auge", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_Auge.png", width=5, height=3)
ind <- ind[-which(ind==150)]
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_Auge_no_outlier.png", width=5, height=3)
autoplot(pca, data = metadata[ind,], colour = 'Group', x=2, y=3)
ggsave("../plots/PCA_Auge_no_outlier_2_3.png", width=5, height=3)

# Create DESeq object
outliers = c(105, 150)
auge = c(141:175)
nonauge = c(1:140)

remove = c(outliers, auge)
remove = c(outliers, nonauge)

remove <- remove[-which(duplicated(remove))]

dds <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~ 0 + Group)
dds <- dds[, -which(dds$SampleNumber %in% remove)]
dds$Group <- factor(metadata$Group[-remove])

vstd <- vst(dds, blind=TRUE)
vstd_mat <- assay(vstd)
colnames(vstd_mat) <- metadata$SampleNumber[-remove]
colnames(vstd_mat)
metadata.short <- metadata[-remove,]

############ Differential analysis

# LRT
design(dds) <- ~ 0 + Group
dds<-DESeq(dds, test = "LRT", full = ~ 0 + Group, reduced = ~1)
resultsNames(dds)

res <- results(dds)
betas <- coef(dds)
colnames(betas)

## LV vs RV vs Sept vs Ap: PEd21

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)


# Wald
dds <- DESeq(dds)
resultsNames(dds)


## LV: PEd21 vs SD21
res <- results(dds, contrast=c("Group","LV_PEd21","LV_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","LV_PEd21","LV_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_PEd21_LV_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_PEd21","LV_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_PEd21_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol
  
pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## LV: PEpp vs SDpp
res <- results(dds, contrast=c("Group","LV_PEpp","LV_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","LV_PEpp","LV_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_PEpp_LV_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_PEpp","LV_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_PEpp_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## LV: np vs SD21
res <- results(dds, contrast=c("Group","LV_np","LV_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","LV_np","LV_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_np_LV_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_np","LV_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_np_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## LV: np vs SDpp
res <- results(dds, contrast=c("Group","LV_np","LV_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","LV_np","LV_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_np_LV_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_np","LV_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_np_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)


## RV: PEd21 vs SD21
res <- results(dds, contrast=c("Group","RV_PEd21","RV_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","RV_PEd21","RV_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/RV_PEd21_RV_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("RV_PEd21","RV_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "RV_PEd21_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## LV: np vs PE21
res <- results(dds, contrast=c("Group","LV_np","LV_PEd21"))
res <- lfcShrink(dds, contrast=c("Group","LV_np","LV_PEd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_np_LV_PEd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_np","LV_PEd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_np_PEd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## LV: np vs PEpp
res <- results(dds, contrast=c("Group","LV_np","LV_PEpp"))
res <- lfcShrink(dds, contrast=c("Group","LV_np","LV_PEpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/LV_np_LV_PEpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("LV_np","LV_PEpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "LV_np_PEpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[-c(1,2),]
resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)



## RV: PEpp vs SDpp
res <- results(dds, contrast=c("Group","RV_PEpp","RV_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","RV_PEpp","RV_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/RV_PEpp_RV_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("RV_PEpp","RV_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "RV_PEpp_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Sept: PEd21 vs SD21
res <- results(dds, contrast=c("Group","Sept_PEd21","Sept_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","Sept_PEd21","Sept_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Sept_PEd21_Sept_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Sept_PEd21","Sept_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Sept_PEd21_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Sept: PEpp vs SDpp
res <- results(dds, contrast=c("Group","Sept_PEpp","Sept_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","Sept_PEpp","Sept_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Sept_PEpp_Sept_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Sept_PEpp","Sept_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Sept_PEpp_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Ap: PEd21 vs SD21
res <- results(dds, contrast=c("Group","Ap_PEd21","Ap_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","Ap_PEd21","Ap_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Ap_PEd21_Ap_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Ap_PEd21","Ap_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Ap_PEd21_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Ap: PEpp vs SDpp
res <- results(dds, contrast=c("Group","Ap_PEpp","Ap_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","Ap_PEpp","Ap_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Ap_PEpp_Ap_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Ap_PEpp","Ap_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]

cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Ap_PEpp_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Auge: PEd21 vs SD21
res <- results(dds, contrast=c("Group","Auge_PEd21","Auge_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","Auge_PEd21","Auge_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_PEd21_Auge_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_PEd21","Auge_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_PEd21_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Auge: PEpp vs SDpp
res <- results(dds, contrast=c("Group","Auge_PEpp","Auge_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","Auge_PEpp","Auge_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_PEpp_Auge_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_PEpp","Auge_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]

cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_PEpp_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Auge: np vs SD21
res <- results(dds, contrast=c("Group","Auge_np","Auge_SDd21"))
res <- lfcShrink(dds, contrast=c("Group","Auge_np","Auge_SDd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_np_Auge_SDd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_np","Auge_SDd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_np_SDd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Auge: np vs SDpp
res <- results(dds, contrast=c("Group","Auge_np","Auge_SDpp"))
res <- lfcShrink(dds, contrast=c("Group","Auge_np","Auge_SDpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_np_Auge_SDpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_np","Auge_SDpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_np_SDpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)


## Auge: np vs PEd21
res <- results(dds, contrast=c("Group","Auge_np","Auge_PEd21"))
res <- lfcShrink(dds, contrast=c("Group","Auge_np","Auge_PEd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_np_Auge_PEd21.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_np","Auge_PEd21"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_np_PEd21"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)

## Auge: np vs SDpp
res <- results(dds, contrast=c("Group","Auge_np","Auge_PEpp"))
res <- lfcShrink(dds, contrast=c("Group","Auge_np","Auge_PEpp"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
write.table(resTable, "../degs/Auge_np_Auge_PEpp.tsv", sep="\t", row.names = F)

ind = which(metadata.short$Group %in% c("Auge_np","Auge_PEpp"))
cluster_vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
cluster_vst_cor <- cor(cluster_vst)
color = hue_pal()(2)
names(color) <- levels(factor(metadata.short[ind,]$Group))
annoCol <- list(Group = color)
txt = "Auge_np_PEpp"
pheatmap(cluster_vst_cor, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_cor_", txt, ".png", sep=""),
         width = 7, height = 5)

resTable <- resTable[which(resTable$symbol!=""),]

vst <- vstd_mat[resTable$ens_gene, as.character(metadata.short$SampleNumber[ind])]
rownames(vst) <- resTable$symbol

pheatmap(vst, annotation = metadata.short[ind,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/WT/heatmap_vst_", txt, ".png", sep=""), scale = "row",
         width = 8, height = 11)


###########################3
## LV vs RV: PEd21
res <- results(dds, contrast=c("Group","LV_PEd21","RV_PEd21"))
res <- lfcShrink(dds, contrast=c("Group","LV_PEd21","RV_PEd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/LV_PEd21_RV_PEd21.tsv", sep="\t", row.names = F)

## RV: PEd21 vs PEpp
res <- results(dds, contrast=c("Group","RV_PEpp","RV_PEd21"))
res <- lfcShrink(dds, contrast=c("Group","RV_PEpp","RV_PEd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/RV_PEpp_RV_PEd21.tsv", sep="\t", row.names = F)

## RV: np vs PEpp
res <- results(dds, contrast=c("Group","RV_PEpp","RV_np"))
res <- lfcShrink(dds, contrast=c("Group","RV_PEpp","RV_np"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/RV_PEpp_RV_np.tsv", sep="\t", row.names = F)

## LV: np vs PEpp
res <- results(dds, contrast=c("Group","LV_PEpp","LV_np"))
res <- lfcShrink(dds, contrast=c("Group","LV_PEpp","LV_np"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/LV_PEpp_LV_np.tsv", sep="\t", row.names = F)

## LV: np vs SDpp
res <- results(dds, contrast=c("Group","LV_SDpp","LV_np"))
res <- lfcShrink(dds, contrast=c("Group","LV_SDpp","LV_np"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/LV_SDpp_LV_np.tsv", sep="\t", row.names = F)

## RV: np vs SDpp
res <- results(dds, contrast=c("Group","RV_SDpp","RV_np"))
res <- lfcShrink(dds, contrast=c("Group","RV_SDpp","RV_np"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/RV_SDpp_RV_np.tsv", sep="\t", row.names = F)

## RV: np vs SDd21
res <- results(dds, contrast=c("Group","RV_SDd21","RV_np"))
res <- lfcShrink(dds, contrast=c("Group","RV_SDd21","RV_np"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]
resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))

write.table(resTable, "../degs/RV_SDd21_RV_np.tsv", sep="\t", row.names = F)


# Heatmap

my_sample_col <- data.frame(ProtocolNum = pm$ProtocolNum)
row.names(my_sample_col) <- colnames(bc)

my_colors <- RColorBrewer::brewer.pal(8, "Set3")
names(my_colors) <- levels(pm$ProtocolNum)
ann_colors = list(ProtocolNum = my_colors)

p <- pheatmap(bc, annotation_col = my_sample_col, annotation_colors = ann_colors)
