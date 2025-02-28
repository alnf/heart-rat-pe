library(CEMiTool)
#library(CeTF)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(viridis)
library(ggrepel)
library(dendextend)
library(ggh4x)
library(dbscan)
library(ggrepel)
library(ggfortify)
library(ggpubr)
library(reshape2)
library(lemon)
library(fgsea)
library("FactoMineR")
library("factoextra")
library(VennDiagram)
source("utils.R")
library(ggVennDiagram)

### Metadata

metadata.norm <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.norm <- metadata.norm[-which(grepl("PE", metadata.norm$Pheno)), ]

### All genes
types = c("_short.tsv", ".tsv")
type = types[1]

region = "LV"

fname = paste(region, "SDd21", region, "SDpp", sep="_")
SDd21_SDpp <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDd21", region, "np", sep="_")
SDd21_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

fname = paste(region, "SDpp", region, "np", sep="_")
SDpp_np <- read.table(paste("../degs/", region, "/degs_", fname, type, sep=""), sep="\t", header = T)

ncolors <- read.table("../metadata/no_colors.tsv", sep="\t", header=T, check.names = F, comment.char = "")
ncols <- ncolors$color
names(ncols) <- ncolors$variable
ncols

### PCA
m <- lcounts[, which(colnames(lcounts) %in% metadata.norm$SampleNumber)]
pca <- prcomp(t(m))
autoplot(pca, data = metadata.norm,
         colour = 'PhenoNames',
         label = TRUE, label.label = "SampleNumber", label.repel=T,
         frame=T, frame.type="norm",
         label.show.legend = F,
         frame.colour = 'PhenoNames') +
  scale_color_manual(values = ncols) +
  scale_fill_manual(values = ncols)
ggsave("../plots/normal_preg/pca.png", width = 8.7, height = 6, scale = 0.8, dpi = 150)
