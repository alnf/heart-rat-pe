library(DESeq2)
library(tximport)
library(yaml)
library(dplyr)
library(biomaRt)
library(ggfortify)
library(RColorBrewer)
library(DEGreport)
library(tibble)
library(pheatmap)
library(scales)

save_pheatmap_png <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

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
metadata$Condition <- sapply(strsplit(metadata$SampleNameUser," "), `[`, 3)
metadata$Condition <- substring(metadata$Condition, 1, 2)

metadata$Pregnancy <- sapply(strsplit(metadata$SampleNameUser," "), `[`, 3)
metadata$Pregnancy <- sub('SD', '', metadata$Pregnancy)
metadata$Pregnancy <- sub('PE', '', metadata$Pregnancy)

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

# Create DESeq object
outliers = c(105, 150)
lv.ind = c(1:35)
rv.ind = c(36:70)
sept.ind = c(71:105)
ap.ind = c(106:140)
auge.ind = c(141:175)

dds.all <- DESeqDataSetFromTximport(txi.kallisto, metadata, ~1)

# LV
txt = "LV"
ind <- lv.ind[!(lv.ind %in% outliers)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
metadata.short$Group <- factor(metadata.short$Group, levels=c("LV_np", "LV_SDd21", "LV_SDpp", "LV_PEd21", "LV_PEpp"))
dds$Group <- metadata.short$Group

# RV
txt = "RV"
ind <- rv.ind[!(rv.ind %in% outliers)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
metadata.short$Group <- factor(metadata.short$Group, levels=c("RV_np", "RV_SDd21", "RV_SDpp", "RV_PEd21", "RV_PEpp"))
dds$Group <- metadata.short$Group

# Sept
txt = "Sept"
ind <- sept.ind[!(sept.ind %in% outliers)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
metadata.short$Group <- factor(metadata.short$Group, levels=c("Sept_np", "Sept_SDd21", "Sept_SDpp", "Sept_PEd21", "Sept_PEpp"))
dds$Group <- metadata.short$Group

# Ap
txt = "Ap"
ind <- ap.ind[!(ap.ind %in% outliers)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
metadata.short$Group <- factor(metadata.short$Group, levels=c("Ap_np", "Ap_SDd21", "Ap_SDpp", "Ap_PEd21", "Ap_PEpp"))
dds$Group <- metadata.short$Group

# Auge
txt = "Auge"
ind <- auge.ind[!(auge.ind %in% outliers)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
metadata.short$Group <- factor(metadata.short$Group, levels=c("Auge_np", "Auge_SDd21", "Auge_SDpp", "Auge_PEd21", "Auge_PEpp"))
dds$Group <- metadata.short$Group

# Heart, np
txt = "heart_np"
ind <- which(metadata$Condition=="np")
ind <- ind[!(ind %in% outliers)]
ind <- ind[!(ind %in% auge.ind)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
#metadata.short$Group <- factor(metadata.short$Group, levels=c("Auge_np", "Auge_SDd21", "Auge_SDpp", "Auge_PEd21", "Auge_PEpp"))
metadata.short$Group <- factor(metadata.short$Group)
dds$Group <- metadata.short$Group

# Heart, pp, SD
txt = "heart_pp_SD"
ind <- which(metadata$Pregnancy=="pp" & metadata$Condition=="SD")
ind <- ind[!(ind %in% outliers)]
ind <- ind[!(ind %in% auge.ind)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
#metadata.short$Group <- factor(metadata.short$Group, levels=c("Auge_np", "Auge_SDd21", "Auge_SDpp", "Auge_PEd21", "Auge_PEpp"))
metadata.short$Group <- factor(metadata.short$Group)
dds$Group <- metadata.short$Group

# Heart, pp, PE
txt = "heart_pp_PE"
ind <- which(metadata$Pregnancy=="pp" & metadata$Condition=="PE")
ind <- ind[!(ind %in% outliers)]
ind <- ind[!(ind %in% auge.ind)]
dds <- dds.all[, which(dds.all$SampleNumber %in% ind)]
metadata.short <- metadata[ind,]
#metadata.short$Group <- factor(metadata.short$Group, levels=c("Auge_np", "Auge_SDd21", "Auge_SDpp", "Auge_PEd21", "Auge_PEpp"))
metadata.short$Group <- factor(metadata.short$Group)
dds$Group <- metadata.short$Group

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
color = hue_pal()(4)
names(color) <- levels(metadata.short$Group)
dev.off()
plotPCA(rld, intgroup="Group") + scale_color_manual(values=color)
ggsave(paste("../plots/PCArlog/PCA_rlog_", txt, ".png", sep=""), width=5, height=3)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
colnames(rld_mat) <- metadata.short$SampleNumber
rld_cor <- cor(rld_mat)

# Plot heatmap
annoCol <- list(Group = color)
pheatmap(rld_cor, annotation = metadata.short[,c("Group"), drop=F], annotation_colors = annoCol,
              filename = paste("../plots/LRT/heatmap_rlog_full_", txt, ".png", sep=""),
              width = 7, height = 5)

# LRT
design(dds) <- ~ 0 + Group
dds<-DESeq(dds, test = "LRT", full = ~ 0 + Group, reduced = ~1)
resultsNames(dds)
res_LRT <- results(dds)

# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.05
padj.cutoff = 0.001
sigLRT_genes <- res_LRT_tb %>% 
  filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)

clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n=1500)

# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
cluster_rlog_cor <- cor(cluster_rlog)
pheatmap(cluster_rlog_cor, annotation = metadata.short[,c("Group"), drop=F], annotation_colors = annoCol,
         filename = paste("../plots/LRT/heatmap_rlog_top1500_", txt, ".png", sep=""),
         width = 7, height = 5)

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = metadata.short, time = "Group", col = NULL, plot = F)
#clusters$plot + geom_boxplot(aes(fill=Condition, color=Condition))
clusters$plot + geom_boxplot(aes(fill=Group, color=Group))

dev.off()
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = metadata.short, time = "Group", col = NULL, plot = F,
                        pattern=c(0,1,0,0))
clusters <- degPatterns(cluster_rlog, metadata = metadata.short, time = "Group", col = NULL, plot = F,
                        pattern=c(1,0,0,0))
clusters <- degPatterns(cluster_rlog, metadata = metadata.short, time = "Group", col = NULL, plot = F,
                        pattern=c(0,0,1,0))
clusters <- degPatterns(cluster_rlog, metadata = metadata.short, time = "Group", col = NULL, plot = F,
                        pattern=c(0,0,0,1))

#clusters$plot + geom_boxplot(aes(fill=Condition, color=Condition))
clusters$plot + geom_boxplot(aes(fill=Group, color=Group))



ggsave(paste("../plots/LRT/clusters_rlog_top1500_", txt, ".png", sep=""), width=16, height=12)

#ggplot(clusters[["normalized"]],
#       aes(Group, value, color = Condition, fill = Condition)) +
#       geom_boxplot() +
#       geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
#       geom_smooth(aes(group=Group), method = "lm")

cluster_groups <- clusters$df
colnames(cluster_groups)[1] <- "ens_gene"
cluster_groups <- merge(cluster_groups, t2g, by="ens_gene")
cluster_groups <- cluster_groups[!duplicated(cluster_groups$ens_gene),]

write.table(cluster_groups[order(cluster_groups$cluster, decreasing = F),],
            paste("../degs/LRT/clusters_rlog_top1500_", txt, ".tsv", sep=""),
            sep="\t", row.names = F)

group4 <- cluster_groups %>%
  filter(cluster == 4)

hex_codes <- hue_pal()(4)
hex_codes 
