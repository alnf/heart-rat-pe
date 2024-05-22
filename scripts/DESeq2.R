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
metadata.h$Region <- factor(sapply(strsplit(metadata.h$Group,"_"), `[`, 1))
metadata.h$Pheno <- factor(sapply(strsplit(metadata.h$Group,"_"), `[`, 2))

dds.h$Group <- factor(metadata.h$Group)

exprs.h.vst <- vst(dds.h, blind=TRUE)
exprs.h <- assay(exprs.h.vst)
colnames(exprs.h) <- metadata.h$SampleNumber
colnames(exprs.h)

### Differential expression in heart
dds.h <- DESeq(dds.h)
resultsNames(dds.h)
counts.dds.h <- counts(dds.h,normalized=TRUE)
lcounts <- counts.dds.h
colnames(lcounts) <- rownames(metadata.h)
lcounts <- log2(lcounts+1)

source("utils.R")

region = "LV"

wrapDegs <- function(dds, t2g, contrast, region, sval.filter=TRUE) {
  contrast[2:3] <- paste(region, contrast[2:3], sep="_")
  fname = paste(contrast[2:3], collapse="_")
  fname_short = paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep="")
  fname = paste("../degs/WT/", region, "/degs_", fname, ".tsv", sep="")
  getDEGs(dds.h, contrast, t2g, sval.filter=sval.filter, filename=fname, filename_short=fname_short)  
}

contrast = c("Group","PEd21","SDd21")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","PEpp","SDpp")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","SDd21","np")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","SDpp","np")
wrapDegs(dds.h, t2g, contrast, region, sval.filter=FALSE)

contrast = c("Group","PEd21","np")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","PEpp","np")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","SDd21","SDpp")
wrapDegs(dds.h, t2g, contrast, region)

contrast = c("Group","PEd21","PEpp")
wrapDegs(dds.h, t2g, contrast, region)

### All genes
types = c("_short.tsv", ".tsv")
type = types[1]

fname = paste(region, "PEd21", region, "SDd21", sep="_")
PEd21_SDd21 <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "SDpp", sep="_")
PEpp_SDpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "PEpp", sep="_")
PEd21_PEpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "SDpp", sep="_")
SDd21_SDpp <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "np", sep="_")
PEd21_np <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "np", sep="_")
PEpp_np <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "np", sep="_")
SDd21_np <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDpp", region, "np", sep="_")
SDpp_np <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

all_genes <- c()
all_genes$ens_gene <- unique(c(PEd21_SDd21$ens_gene, PEpp_SDpp$ens_gene, PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene,
                               PEd21_np$ens_gene, PEpp_np$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene))
all_genes <- as.data.frame(all_genes)

### Heatmaps

#TODO
source("utils.R")
contrast = c("SDd21", "PEd21")
fname <- paste(region, contrast, sep="_")
fname = paste(fname, collapse="_")
fname = paste("../plots/WT/", region, "/heatmap_vst_", fname, ".png", sep="")
makeHeatmap(dds.h, contrast, PEd21_SDd21, lFC=2, fname)

### Venn diagram: SD vs PE
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/WT/", region, "/venn_", fname, ".png", sep="")
makeVenn(2, list(PEd21_SDd21$symbol, PEpp_SDpp$symbol), names, region,
         fname, brewer.pal(3, "Pastel2")[1:2])

dd <- intersect(PEd21_SDd21$symbol, PEpp_SDpp$symbol)
ss <- PEd21_SDd21$symbol[!(PEd21_SDd21$symbol %in% dd)]

### Single gene plots
source("utils.R")
region = "LV"
groups <- dds.h$Group[grepl(region, dds.h$Group)]
gene = ""
title = unique(t2g$symbol[which(t2g$ens_gene==gene)])
my_comparisons <- list(c("LV_PEd21", "LV_SDd21"), c("LV_PEpp", "LV_SDpp"))
fname = paste("../plots/WT/", region, "/gene_", title, ".png", sep="")
genePlot(dds.h, gene, "Group", groups, my_comparisons, fname, title)

### Region heatmap
source("utils.R")
filename = paste("../plots/WT/", region, "/heatmap_vst_cor_", region, ".png", sep="")
regionPlot(dds.h, region, all_genes$ens_gene, filename)

### Venn between regions
comparisons <- c("PEd21", "SDd21")

region = "LV"
fname = paste(region, comparisons[1], region, comparisons[2], sep="_")
LV <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

region = "RV"
fname = paste(region, comparisons[1], region, comparisons[2], sep="_")
RV <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

region = "Sept"
fname = paste(region, comparisons[1], region, comparisons[2], sep="_")
Sept <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

region = "Ap"
fname = paste(region, comparisons[1], region, comparisons[2], sep="_")
Ap <- read.table(paste("../degs/WT/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

names = c("LV", "RV", "Sept", "Ap")
fname = paste(comparisons, collapse="_")
all_genes <- unique(c(LV$symbol, RV$symbol, Sept$symbol, Ap$symbol))
all_genes <- unique(c(LV$ens_gene, RV$ens_gene, Sept$ens_gene, Ap$ens_gene))
#all_genes <- unique(c(LV$symbol, Sept$symbol, Ap$symbol))

title = paste(fname, ", ngenes=", length(all_genes), sep="")
makeVenn(4, list(LV$symbol, RV$symbol, Sept$symbol, Ap$symbol), names, title,
         paste("../plots/WT/venn_regions_", fname, ".png", sep=""),
         brewer.pal(4, "Pastel2"), dist = 0.1)

#makeVenn(3, list(LV$symbol, Sept$symbol, Ap$symbol), c("LV", "Sept", "Ap"), paste("LV:", title),
#         paste("../plots/WT/venn_LV_", fname, ".png", sep=""),
#         brewer.pal(3, "Pastel2"), dist = 0.1)

### Heatmap of intersection genes

source("utils.R")
int_genes <- c()
int_genes$ens_gene <- Reduce(intersect, list(LV$ens_gene,
                                             RV$ens_gene,
                                             Sept$ens_gene, Ap$ens_gene))
int_genes$symbol <- LV$symbol[which(LV$ens_gene %in% int_genes$ens_gene)]
int_genes <- as.data.frame(int_genes)
rownames(int_genes) <- int_genes$ens_gene
  
LV <- LV[LV$ens_gene %in% int_genes$ens_gene, c("ens_gene", "log2FoldChange")]
RV <- RV[RV$ens_gene %in% int_genes$ens_gene, c("ens_gene", "log2FoldChange")]
Sept <- Sept[Sept$ens_gene %in% int_genes$ens_gene, c("ens_gene", "log2FoldChange")]
Ap <- Ap[Ap$ens_gene %in% int_genes$ens_gene, c("ens_gene", "log2FoldChange")]

df_list <- list(LV,
                RV,
                Sept, Ap)
merged <- df_list %>% reduce(full_join, by='ens_gene')
rownames(merged) <- merged$ens_gene
merged <- merged[,2:ncol(merged)]
colnames(merged) <- c("LV",
                      "RV",
                      "Sept", "Ap")

write.table(int_genes, "../degs/common_between_3_regions_pp.tsv", sep="\t", row.names = F, col.names = T)


ind = which(rowSds(as.matrix(merged)) > 0.3)
merged <- merged[ind,]

fname = paste(comparisons, collapse = "_")
#filename = paste("../plots/WT/heatmap_vst_cor_all_", fname, ".png", sep="")
filename = paste("../plots/WT/heatmap_vst_cor_LV_", fname, ".png", sep="")

makeHeatmapFC(merged, int_genes[rownames(merged),], comparisons, filename, "none")
makeHeatmapFC(merged, int_genes[rownames(merged),], comparisons, filename, "row")

int_genes_a <- int_genes
int_genes_b <- int_genes

common <- intersect(int_genes_a$ens_gene, int_genes_b$ens_gene)
common <- int_genes_a[common,]
write.table(common, "../degs/common_between_all_regions_condition.tsv", sep="\t", row.names = F, col.names=T)
#write.table(common, "../degs/common_between_3_regions_condition.tsv", sep="\t", row.names = F, col.names = T)
