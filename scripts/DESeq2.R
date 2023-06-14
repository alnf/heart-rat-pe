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

region = "Ap"
contrast = c("Group","PEd21","SDd21")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","PEpp","SDpp")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","SDd21","np")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, lFCvis = 1, width=7, height=5)

contrast = c("Group","SDpp","np")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
#For RV and Ap
#getDEGs(dds.h, contrast, t2g, lFCvis = 1, width=7, height=5, sval.filter=FALSE)
getDEGs(dds.h, contrast, t2g, lFCvis = 1, width=7, height=5)

contrast = c("Group","PEd21","np")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","PEpp","np")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","SDd21","SDpp")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

contrast = c("Group","PEd21","PEpp")
contrast[2:3] <- paste(region, contrast[2:3], sep="_")
getDEGs(dds.h, contrast, t2g, width=7, height=5)

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

# Venn diagram: SD vs PE
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/WT/", region, "/venn_", fname, ".png", sep="")
makeVenn(2, list(PEd21_SDd21$symbol, PEpp_SDpp$symbol), names, region,
         fname, brewer.pal(3, "Pastel2")[1:2])

dd <- intersect(PEd21_SDd21$symbol, PEpp_SDpp$symbol)
ss <- PEd21_SDd21$symbol[!(PEd21_SDd21$symbol %in% dd)]

### Single gene plots
source("utils.R")
groups <- dds.h$Group[grepl(region, dds.h$Group)]
#gene = "ENSRNOG00000048222"
gene = "ENSRNOG00000016451"
title = unique(t2g$symbol[which(t2g$ens_gene==gene)])
my_comparisons <- list(c("LV_PEd21", "LV_SDd21"), c("LV_PEpp", "LV_SDpp"))
fname = paste("../plots/WT/", region, "/gene_", title, ".png", sep="")
genePlot(dds.h, gene, "Group", groups, my_comparisons, fname, title)

### Region heatmap
source("utils.R")
regionPlot(dds.h, region, all_genes$ens_gene)
#regionPlot(dds.h, region, rownames(dds.h))


### Venn between regions
comparisons <- c("PEd21", "SDd21")
comparisons <- c("PEpp", "SDpp")

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
title = paste(fname, ", ngenes=", length(all_genes), sep="")
makeVenn(4, list(LV$symbol, RV$symbol, Sept$symbol, Ap$symbol), names, title,
         paste("../plots/WT/venn_regions_", fname, ".png", sep=""),
         brewer.pal(4, "Pastel2"), dist = 0.1)



