library(DESeq2)
library(tximport)
library(yaml)
library(dplyr)
library(biomaRt)
library(ggfortify)
library(RColorBrewer)
library(VennDiagram)

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
                                     "rgd_symbol", "mgi_symbol", "entrezgene_description"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, symbol = rgd_symbol, msymbol = mgi_symbol, description=entrezgene_description)

# Import kallisto data
files <- file.path(kal_dirs, "abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", tx2gene=t2g[,1:2], ignoreTxVersion = T, txOut = F)
counts <- txi.kallisto$counts
head(counts)

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
outliers = c(7, 105, 150)
auge = c(141:175)
nonauge = c(1:140)

auge.ind = c(outliers, auge)
auge.ind <- auge.ind[-which(duplicated(auge.ind))]
heart.ind = c(outliers, nonauge)
heart.ind <- heart.ind[-which(duplicated(heart.ind))]

dds <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~ 0 + Group)
#dds.back <- dds

######### Heart analysis
dds.h <- dds[, -which(dds$SampleNumber %in% auge.ind)]
metadata.h <- metadata[-auge.ind,]
dds.h$Group <- factor(metadata.h$Group)

exprs.h.vst <- vst(dds.h, blind=TRUE)
exprs.h <- assay(exprs.h.vst)
colnames(exprs.h) <- metadata.h$SampleNumber
colnames(exprs.h)

### Differential expression in heart
dds.h <- DESeq(dds.h)
resultsNames(dds.h)
counts.dds.h <- counts(dds.h,normalized=TRUE)

source("utils.R")

contrast = c("Group","LV_PEd21","LV_SDd21")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","LV_PEpp","LV_SDpp")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","LV_SDd21","LV_np")
getDEGs(dds.h, contrast, t2g, lFCvis = 1, width=7, height=5)

contrast = c("Group","LV_SDpp","LV_np")
getDEGs(dds.h, contrast, t2g, lFCvis = 1, width=7, height=5)

contrast = c("Group","LV_PEd21","LV_np")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","LV_PEpp","LV_np")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","LV_SDd21","LV_SDpp")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","LV_PEd21","LV_PEpp")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

### Single gene plots

my_comparisons <- list(c("preterm34weeksnoFGR", "term37weeksnoFGR"))
ggboxplot(data=pdata.ol, x="Condition", y="bmi", color="Condition", palette = "jco", add = "jitter") +
  stat_compare_means(method = "t.test", comparisons = my_comparisons)



### Venn diagrams
region = "LV"

fname = paste(region, "PEd21", region, "SDd21", sep="_")
PEd21_SDd21 <- read.table(paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "SDpp", sep="_")
PEpp_SDpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "PEpp", sep="_")
PEd21_PEpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "SDpp", sep="_")
SDd21_SDpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep=""), sep="\t", header = T)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(PEd21_SDd21$symbol, PEpp_SDpp$symbol),
  category.names = c("PEd21_SDd21" , "PEpp_SDpp"),
  filename = '../plots/WT/LV/venn_diagramm1.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-5, 5),
  cat.dist = c(0.055, 0.055)
)

myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(SDd21_SDpp$symbol, PEd21_PEpp$symbol),
  category.names = c("SDd21_SDpp" , "PEd21_PEpp"),
  filename = '../plots/WT/LV/venn_diagramm2.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[3:4],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-5, 5),
  cat.dist = c(0.055, 0.055)
)

dd <- intersect(SDd21_SDpp$symbol, PEd21_PEpp$symbol)


dd <- intersect(PEd21_SDd21$symbol, PEpp_SDpp$symbol)

ss <- PEpp_SDpp$symbol[!(PEpp_SDpp$symbol %in% dd)]
######### Eye analysis




#########
dds <- dds[, -which(dds$SampleNumber %in% remove)]
dds$Group <- factor(metadata$Group[-remove])

vstd <- vst(dds, blind=TRUE)
vstd_mat <- assay(vstd)
colnames(vstd_mat) <- metadata$SampleNumber[-remove]
colnames(vstd_mat)
metadata.short <- metadata[-remove,]

############ Differential analysis

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
