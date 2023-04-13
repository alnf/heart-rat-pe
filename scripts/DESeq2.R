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
dds <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~Group)

# Differential analysis
dds$Group <- relevel(dds$Group, ref = "LV_np")
dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients
resultsNames(dds)[grepl('Group_LV', resultsNames(dds))]
res <- results(dds, name="Group_LV_SDd21_vs_LV_np")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="Group_LV_SDd21_vs_LV_np", type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]

###
resultsNames(dds)[grepl('Group_LV', resultsNames(dds))]
res <- results(dds, name="Group_LV_PEd21_vs_LV_np")
res <- lfcShrink(dds, coef="Group_LV_PEd21_vs_LV_np", type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 2)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]


###
res <- results(dds, contrast=c("Group","LV_SDd21","LV_PEd21"))
res <- lfcShrink(dds, contrast=c("Group","LV_SDd21","LV_PEd21"), type="ashr")

resTable <- data.frame(res)
resTable <- resTable %>% filter(padj <= 0.05)
resTable <- resTable %>% filter(abs(log2FoldChange) > 1)

resTable$ens_gene <- rownames(resTable)

resTable <- merge(resTable, t2g, by="ens_gene")
resTable <- resTable[!duplicated(resTable$ens_gene),]

