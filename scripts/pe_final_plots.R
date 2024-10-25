### Metadata

metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$Pheno == "np"), ]

### All genes
types = c("_short.tsv", ".tsv")
type = types[1]

region = "LV"
fname = paste(region, "PEd21", region, "SDd21", sep="_")
PEd21_SDd21 <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEpp", region, "SDpp", sep="_")
PEpp_SDpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "PEd21", region, "PEpp", sep="_")
PEd21_PEpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "SDpp", sep="_")
SDd21_SDpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

PEd21_SDd21 <- PEd21_SDd21[-which(PEd21_SDd21$ens_gene==agt$ens_gene),]
PEpp_SDpp <- PEpp_SDpp[-which(PEpp_SDpp$ens_gene==agt$ens_gene),]

mcolors <- read.table("../metadata/pe_colors.tsv", sep="\t", header=T, check.names = F, comment.char = "")
mcols <- mcolors$color
names(mcols) <- mcolors$variable

### PCA
m <- lcounts[, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
pca <- prcomp(t(m))
autoplot(pca, data = metadata.left,
         colour = 'PhenoNames',
         label = TRUE, label.label = "SampleNumber", label.repel=T,
         frame=T, frame.type="norm",
         label.show.legend = F,
         frame.colour = 'PhenoNames') +
       scale_color_manual(values = mcols) +
       scale_fill_manual(values = mcols)
ggsave("../plots/final/pca.png", width = 8.7, height = 6, scale = 0.8, dpi = 150)


## Plotting Venns

### preg vs post
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/final/final_venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost")
makeVenn(2, list(PEd21_SDd21$ens_gene[which(abs(PEd21_SDd21$log2FoldChange)>1)],
                 PEpp_SDpp$ens_gene[which(abs(PEpp_SDpp$log2FoldChange)>1)]),
         names,
         paste("DEGs in pregnancy and postpartum between PE and WT,\n lFC=1, ngenes=", nrow(genes_pe_reduced_lfc1), sep=""),
         fname, c(mcols["PEpregWT"], mcols["PEpostWT"]), width=650, height=600)


### preg, post, preg_effect
names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/final/final_venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost", "PEpreg_PEpost", "WTpreg_WTpost")
makeVenn(4, list(PEd21_SDd21$ens_gene[which(abs(PEd21_SDd21$log2FoldChange)>1)], PEpp_SDpp$ens_gene[which(abs(PEpp_SDpp$log2FoldChange)>1)],
                 PEd21_PEpp$ens_gene[which(abs(PEd21_PEpp$log2FoldChange)>1)], SDd21_SDpp$ens_gene[which(abs(SDd21_SDpp$log2FoldChange)>1)]),
         names,
         paste("DEGs in pregnancy and postpartum between PE and WT,\n and pregnancy in both PE and WT,\n lFC=1, ngenes=", nrow(genes_pe_lfc1), sep=""),
         fname, c(mcols["PEpregWT"], mcols["PEpostWT"], mcols["pregPEpost"], mcols["pregWTpost"]), dist=0.15, width=850, height=700)


## Plotting heatmaps, ora

### preg vs post

### Calculate modules
anno_data <- metadata.left[,c("SampleNumber", "PhenoNames")]
anno_data$SampleNumber <- as.character(anno_data$SampleNumber)
colnames(anno_data) <- c("SampleName", "Class")
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
m <- lcounts[genes_pe_reduced_lfc1$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

cem <- cemitool(as.data.frame(m), anno_data, filter = F,
                force_beta = T, apply_vst = F, rank_method = "median",
                gsea_min_size = 10
)
cem
cem@enrichment_plot

### Plot heatmap
m1 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M1"),]$genes)
m2 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M2"),]$genes)
m3 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M3"),]$genes)

dlists <- list(preg=PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),], post=PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),])
cols <- list(preg=mcols["PEpregWT"], post=mcols["PEpostWT"])
acols <- list(Group = mcols)
modcols <- mcols[c((9:11),14)]

rspl <- data.frame(modules = cem@module$modules)
rspl$ens_gene <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M3", "NC"))

glist <- genes_pe_reduced_lfc1
glist$ens_gene <- genes_pe_reduced_lfc1$ens_gene[order(match(genes_pe_reduced_lfc1$ens_gene, rspl$ens_gene))]

# Without clusteting
order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=F, csplit=T, dlists, cols, acols, clegend=F, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_reduced_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()

# With clusterting
order <- c("PEpreg", "PEpost", "WTpreg", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=T, csplit=T, dlists, cols, acols, clegend=T, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_reduced_lfc1_dd.png", width = 10.5, height = 9.1, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()

### preg, post, preg_effect

### Calculate modules
anno_data <- metadata.left[,c("SampleNumber", "PhenoNames")]
anno_data$SampleNumber <- as.character(anno_data$SampleNumber)
colnames(anno_data) <- c("SampleName", "Class")
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
m <- lcounts[genes_pe_lfc1$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

cem <- cemitool(as.data.frame(m), anno_data, filter = F,
                force_beta = T, apply_vst = F, rank_method = "median",
                gsea_min_size = 10
)
cem
cem@enrichment_plot

### Plot heatmap

m1 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M1"),]$genes)
m2 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M2"),]$genes)
m3 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M3"),]$genes)
m4 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M4"),]$genes)
m5 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M5"),]$genes)

dlists <- list(preg=PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),], post=PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),],
               pregPEpost=PEd21_PEpp[which(abs(PEd21_PEpp$log2FoldChange)>1),], pregWTpost=SDd21_SDpp[which(abs(SDd21_SDpp$log2FoldChange)>1),])
cols <- list(preg=mcols["PEpregWT"], post=mcols["PEpostWT"], pregPEpost=mcols["pregPEpost"], pregWTpost=mcols["pregWTpost"])
acols <- list(Group = mcols)
modcols <- mcols[9:14]

rspl <- data.frame(modules = cem@module$modules)
rspl$ens_gene <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M3", "M4", "M5","NC"))

glist <- genes_pe_lfc1
glist$ens_gene <- genes_pe_lfc1$ens_gene[order(match(genes_pe_lfc1$ens_gene, rspl$ens_gene))]

# Without clusteting
order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=F, csplit=T, dlists, cols, acols, clegend=F, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()

# With clusterting
order <- c("PEpreg", "PEpost", "WTpreg", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=T, csplit=T, dlists, cols, acols, clegend=T, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_lfc1_dd.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()
