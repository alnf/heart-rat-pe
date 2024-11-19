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
gene = "ENSRNOG00000013578"

title = unique(t2g$symbol[which(t2g$ens_gene==gene)])
my_comparisons <- list(c("LV_PEd21", "LV_SDd21"), c("LV_PEpp", "LV_SDpp"))
fname = paste("../plots/genes/gene_", title, ".png", sep="")
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
