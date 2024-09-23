library(DESeq2)
library(tximport)
library(yaml)
library(dplyr)
library(biomaRt)
library(ggfortify)
library(RColorBrewer)
library(VennDiagram)
library(tidyverse)
source("utils.R")

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
metadata$GroupNames <- metadata$Group
metadata$GroupNames <- gsub("d21", "preg", metadata$GroupNames)
metadata$GroupNames <- gsub("np", "SDnp", metadata$GroupNames)
metadata$GroupNames <- gsub("SD", "WT", metadata$GroupNames)
metadata$GroupNames <- gsub("pp", "post", metadata$GroupNames)

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
ggsave("../plots/QC/PCA_all.png", width=7, height=4, dpi=200)

ind <- grep("Auge", metadata$Group, invert=T)
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart.png", width=7, height=5, dpi=180)

ind <- which(!grepl("Auge", metadata$Group) & grepl("np", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_np.png", width=5, height=3, dpi=200)

ind <- which(!grepl("Auge", metadata$Group) & grepl("pp", metadata$Group) & grepl("SD", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_SDpp.png", width=5, height=3, dpi=200)

ind <- which(!grepl("Auge", metadata$Group) & grepl("d21", metadata$Group) & grepl("SD", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_SDd21.png", width=5, height=3, dpi=200)

ind <- which(!grepl("Auge", metadata$Group) & grepl("pp", metadata$Group) & grepl("PE", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_PEpp.png", width=5, height=3, dpi=200)

# Remove outlier
ind <- ind[-which(ind==105)]
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_PEpp_no_outlier.png", width=5, height=3, dpi=200)

ind <- which(!grepl("Auge", metadata$Group) & grepl("d21", metadata$Group) & grepl("PE", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_heart_PEd21.png", width=5, height=3, dpi=200)

ind <- which(grepl("LV", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_LV.png", width=5, height=3, dpi=200)

ind <- which(grepl("RV", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_RV.png", width=5, height=3, dpi=200)

ind <- which(grepl("Ap", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_Ap.png", width=5, height=3, dpi=200)

ind <- which(grepl("Sept", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_Sept.png", width=5, height=3, dpi=200)

ind <- which(grepl("Auge", metadata$Group))
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/PCA_Auge.png", width=5, height=3, dpi=200)

ind <- ind[-which(ind==150)]
pca <- prcomp(t(log2(txi.kallisto$counts[,ind]+0.5)))
autoplot(pca, data = metadata[ind,], colour = 'Group')
ggsave("../plots/QC/PCA_Auge_no_outlier.png", width=5, height=3, dpi=200)
autoplot(pca, data = metadata[ind,], colour = 'Group', x=2, y=3)
ggsave("../plots/QC/PCA_Auge_no_outlier_PC23.png", width=5, height=3, dpi=200)

# Create DESeq object
outliers = c(7, 105, 150)
#outliers = c(7, 105, 150, 13, 18)
auge = c(141:175)
heart = c(1:140)

auge.ind = which(!(auge %in% outliers))
heart.ind <- which(!(heart %in% outliers))

dds <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~ 0 + Group)

######### Heart analysis
dds.h <- dds[, which(dds$SampleNumber %in% heart.ind)]
metadata.h <- metadata[heart.ind,]
metadata.h$Region <- factor(sapply(strsplit(metadata.h$Group,"_"), `[`, 1))
metadata.h$Pheno <- factor(sapply(strsplit(metadata.h$Group,"_"), `[`, 2))
metadata.h$PhenoNames <- factor(sapply(strsplit(metadata.h$GroupNames,"_"), `[`, 2))


dds.h$Group <- factor(metadata.h$Group)

exprs.h.vst <- vst(dds.h, blind=TRUE)
exprs.h <- assay(exprs.h.vst)
colnames(exprs.h) <- metadata.h$SampleNumber
colnames(exprs.h)

### Differential expression in heart
dds.h <- DESeq(dds.h)
# sizeFactors is supposed to be NULL when using tximport https://support.bioconductor.org/p/117781/
normalizationFactors(dds.h)

resultsNames(dds.h)
lcounts <- counts(dds.h,normalized=TRUE)
colnames(lcounts) <- rownames(metadata.h)
lcounts <- log2(lcounts+1)
rcounts <- counts(dds.h,normalized=FALSE)
colnames(rcounts) <- rownames(metadata.h)
tpm_abd <- txi.kallisto$abundance[,heart.ind]
colnames(tpm_abd) <- rownames(metadata.h)

# Writing expression data in files
write.table(lcounts, "../exprs/lcounts_heart.tsv", sep="\t", row.names = T, col.names=T)
write.table(rcounts, "../exprs/rcounts_heart.tsv", sep="\t", row.names = T, col.names=T)
write.table(exprs.h, "../exprs/vst_heart.tsv", sep="\t", row.names = T, col.names=T)
write.table(tpm_abd, "../exprs/tpm_abd_heart.tsv", sep="\t", row.names = T, col.names=T)
save(lcounts, file = "../exprs/lcounts_heart.Rdata")
save(rcounts, file = "../exprs/rcounts_heart.Rdata")
save(exprs.h, file = "../exprs/vst_heart.Rdata")
save(tpm_abd, file = "../exprs/tpm_abd_heart.Rdata")

wrapDegs <- function(dds, t2g, contrast, region, sval.filter=TRUE, lFC=1) {
  contrast[2:3] <- paste(region, contrast[2:3], sep="_")
  fname = paste(contrast[2:3], collapse="_")
  fname_short = paste("../degs/", region, "/degs_", fname, "_short.tsv", sep="")
  fname = paste("../degs/", region, "/degs_", fname, ".tsv", sep="")
  getDEGs(dds.h, contrast, t2g, sval.filter=sval.filter, filename=fname, filename_short=fname_short, lFC=lFC)  
}

region = "LV"

contrast = c("Group","PEd21","SDd21")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","PEpp","SDpp")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","SDd21","np")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","SDpp","np")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","PEd21","np")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","PEpp","np")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","SDd21","SDpp")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

contrast = c("Group","PEd21","PEpp")
wrapDegs(dds.h, t2g, contrast, region, lFC=0.58)

### All genes
types = c("_short.tsv", ".tsv")
type = types[1]

fname = paste(region, "PEd21", region, "SDd21", sep="_")
PEd21_SDd21 <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "SDpp", sep="_")
PEpp_SDpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "PEpp", sep="_")
PEd21_PEpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "SDpp", sep="_")
SDd21_SDpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "np", sep="_")
PEd21_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "np", sep="_")
PEpp_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "np", sep="_")
SDd21_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDpp", region, "np", sep="_")
SDpp_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

# Getting unions of genes
agt <- t2g[which(t2g$symbol == "Agt"),]

genes_pe <- c()
genes_pe$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene))
genes_pe <- as.data.frame(genes_pe)
genes_pe <- genes_pe[-which(genes_pe$ens_gene==agt$ens_gene),,drop=F]

genes_pe_reduced <- c()
genes_pe_reduced$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene))
genes_pe_reduced <- as.data.frame(genes_pe_reduced)
genes_pe_reduced <- genes_pe_reduced[-which(genes_pe_reduced$ens_gene==agt$ens_gene),,drop=F]

genes_norm <- c()
genes_norm$ens_gene <- unique(c(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene))
genes_norm <- as.data.frame(genes_norm)

## Plotting venns and heatmaps

### preg vs post
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost")
makeVenn(2, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)]),
         names,
         paste("DEGs in pregnancy and postpartum,\n lFC=0.58, ngenes=", nrow(genes_pe_reduced), " (genes_pe_reduced list)", sep=""),
         fname, brewer.pal(3, "Pastel2")[1:2], width=600, height=550)

pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp)
cols <- list(preg="pink", post="lightblue")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_reduced.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=T, dlists, cols, annoCol, clegend=F)
png("../plots/heatmap_genes_pe_reduced_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_reduced_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

### preg, post, preg_effect
names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost", "PEpreg_PEpost", "WTpreg_WTpost")
makeVenn(4, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)],
                 PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene),
         names,
         paste("DEGs in pregnancy and postpartum and pregnancy effect,\n lFC=0.58, ngenes=", nrow(genes_pe), " (genes_pe list)", sep=""),
         fname, brewer.pal(4, "Pastel2"), dist=0.15, width=800, height=650)


pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp, pe_pref=PEd21_PEpp, wt_pref=SDd21_SDpp)
cols <- list(preg="pink", post="lightblue", pe_pref="orange", wt_pref="maroon")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=T, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=T, dlists, cols, annoCol, clegend=T, clcol=F)
png("../plots/heatmap_genes_pe_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

### norm pregnancy
names = c("SDd21_SDpp", "SDd21_np", "SDpp_np")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, ".png", sep="")
names = c("WTpreg_WTpost", "WTpreg_WTnp", "WTpost_WTnp")
makeVenn(3, list(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene), names,
         paste("DEGs in normal pregnancy,\n lFC=0.58, ngenes=", nrow(genes_norm), " (genes_norm list)", sep=""),
         fname, brewer.pal(3, "Pastel2"))

## lFC = 1
PEd21_SDd21 <- PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),]
PEpp_SDpp <- PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),]
PEd21_PEpp <- PEd21_PEpp[which(abs(PEd21_PEpp$log2FoldChange)>1),]
SDd21_SDpp <- SDd21_SDpp[which(abs(SDd21_SDpp$log2FoldChange)>1),]
SDd21_np <- SDd21_np[which(abs(SDd21_np$log2FoldChange)>1),]
SDpp_np <- SDpp_np[which(abs(SDpp_np$log2FoldChange)>1),]

genes_pe_lfc1 <- c()
genes_pe_lfc1$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene))
genes_pe_lfc1 <- as.data.frame(genes_pe_lfc1)
genes_pe_lfc1 <- genes_pe_lfc1[-which(genes_pe_lfc1$ens_gene==agt$ens_gene),,drop=F]

genes_pe_reduced_lfc1 <- c()
genes_pe_reduced_lfc1$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene))
genes_pe_reduced_lfc1 <- as.data.frame(genes_pe_reduced_lfc1)
genes_pe_reduced_lfc1 <- genes_pe_reduced_lfc1[-which(genes_pe_reduced_lfc1$ens_gene==agt$ens_gene),,drop=F]

genes_norm_lfc1 <- c()
genes_norm_lfc1$ens_gene <- unique(c(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene))
genes_norm_lfc1 <- as.data.frame(genes_norm_lfc1)

## Plotting venns for lFC=1
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost")
makeVenn(2, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)]),
         names,
         paste("DEGs in pregnancy and postpartum,\n lFC=1, ngenes=", nrow(genes_pe_reduced_lfc1), " (genes_pe_reduced_lfc1 list)", sep=""),
         fname, brewer.pal(3, "Pastel2")[1:2], width=600, height=550)

pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp)
cols <- list(preg="pink", post="lightblue")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_reduced_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=T, dlists, cols, annoCol, clegend=F)
png("../plots/heatmap_genes_pe_reduced_lfc1_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_reduced_lfc1_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

###
names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost", "PEpreg_PEpost", "WTpreg_WTpost")
makeVenn(4, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)],
                 PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene),
         names,
         paste("DEGs in pregnancy and postpartum and pregnancy effect,\n lFC=1, ngenes=", nrow(genes_pe_lfc1), " (genes_pe_lfc1 list)", sep=""),
         fname, brewer.pal(4, "Pastel2"), dist=0.15, width=800, height=650)

pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp, pe_pref=PEd21_PEpp, wt_pref=SDd21_SDpp)
cols <- list(preg="pink", post="lightblue", pe_pref="orange", wt_pref="maroon")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=T, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_lfc1_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=T, dlists, cols, annoCol, clegend=T, clcol=F)
png("../plots/heatmap_genes_pe_lfc1_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()


###
names = c("SDd21_SDpp", "SDd21_np", "SDpp_np")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("WTpreg_WTpost", "WTpreg_WTnp", "WTpost_WTnp")
makeVenn(3, list(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene), names,
         paste("DEGs in normal pregnancy,\n lFC=1, ngenes=", nrow(genes_norm_lfc1), " (genes_norm_lfc1 list)", sep=""),
         fname, brewer.pal(3, "Pastel2"), width=800, height=750)


########## Old
all_genes <- c()
all_genes$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene,
                               PEd21_np$ens_gene, PEpp_np$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene))
all_genes <- as.data.frame(all_genes)


all_genes_final_strict <- c()
all_genes_final_strict$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene))
all_genes_final_strict <- as.data.frame(all_genes_final_strict)

PEd21_SDd21 <- PEd21_SDd21[which(PEd21_SDd21$padj<0.05),]
PEpp_SDpp <- PEpp_SDpp[which(PEpp_SDpp$padj<0.05),]
PEd21_PEpp <- PEd21_PEpp[which(PEd21_PEpp$padj<0.05),]
SDd21_SDpp <- SDd21_SDpp[which(SDd21_SDpp$padj<0.05),]

PEd21_SDd21 <- PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),]
PEpp_SDpp <- PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),]
PEd21_PEpp <- PEd21_PEpp[which(abs(PEd21_PEpp$log2FoldChange)>1),]
SDd21_SDpp <- SDd21_SDpp[which(abs(SDd21_SDpp$log2FoldChange)>1),]

all_genes_final <- c()
all_genes_final$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene))
all_genes_final <- as.data.frame(all_genes_final)

all_genes_WT_PE <- c()
all_genes_WT_PE$ens_gene <- unique(c(PEpp_SDpp$ens_gene, PEd21_SDd21$ens_gene))
all_genes_WT_PE <- as.data.frame(all_genes_WT_PE)







