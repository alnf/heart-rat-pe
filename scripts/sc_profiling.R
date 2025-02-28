library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(UCell)

#devtools::install_github('immunogenomics/presto')


########## Huebner heart atlass

so <- readRDS("../metadata/huebner_heart_cell_atlas.rds")
#remove(so)
gc()
someta <- so@meta.data
head(Features(so))
ftrs <- so@assays$RNA@meta.features

sosub <- subset(so, subset = Region_x == "LV")
sosub <- subset(sosub, subset = sex == "female")
sosub <- subset(sosub, subset = tissue == "heart left ventricle")
someta <- sosub@meta.data

sonorm <- subset(sosub, subset = disease == "normal")
sodc <- subset(sosub, subset = disease == "dilated cardiomyopathy")

Idents(sonorm) <- "cell_type"
Idents(sodc) <- "cell_type"


module = "m1"
modules <- c("m1", "m2", "m3", "m4", "m5")
markers <- list(m1=c(), m2=c(), m3=c(), m4=c(), m5=c())

for (i in 1:length(modules)) {
  mlist <- merge(eval(parse(text = modules[i])), t2gh, by="ens_gene")
  mlist <- mlist[-which(mlist$hsapiens_homolog_ensembl_gene==""), ]
  if (length(which(duplicated(mlist$hsapiens_homolog_ensembl_gene)))>0) { 
    mids <- mlist$hsapiens_homolog_ensembl_gene[-which(duplicated(mlist$hsapiens_homolog_ensembl_gene))]
  } else {
    mids <- mlist$hsapiens_homolog_ensembl_gene
  }
    
  markers[[i]] <- mids
}

sonorm <- AddModuleScore_UCell(sonorm, features = markers)
signature.names <- paste0(names(markers), "_UCell")
VlnPlot(sonorm, features = signature.names, group.by = "cell_type", pt.size=0)

sodc <- AddModuleScore_UCell(sodc, features = markers)
signature.names <- paste0(names(markers), "_UCell")
VlnPlot(sodc, features = signature.names, group.by = "cell_type", pt.size=0)




options(future.globals.maxSize= 4891289600)
myeloid <- FindMarkers(sonorm, ident.1 = "myeloid cell", verbose = FALSE) %>%
  arrange(-avg_log2FC) %>%
  rownames_to_column(var = "gene") %>%
  pull(gene) %>% 
  .[1:100]
sonorm <- AddModuleScore(sonorm, features = list(myeloid), name="myeloid")

sonorm <- AddModuleScore(sonorm, features = list(mids), name=module)
sonormmeta <- sonorm@meta.data
FeaturePlot(sonorm, features = paste(module, "1", sep=""), label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  labs(title = paste(toupper(module), "enrichment"))
FeaturePlot(sonorm, features = "myeloid1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  labs(title = paste(toupper(module), "enrichment"))

p <- ggplot(sonormmeta, aes(x=cell_type, y=m11, fill=cell_type)) + 
  geom_violin() +  geom_boxplot(width=0.1)
p <- ggplot(sonormmeta, aes(x=cell_type, y=myeloid1, fill=cell_type)) + 
  geom_violin() +  geom_boxplot(width=0.1)
p <- ggplot(sonormmeta, aes(x=cell_type, y=m5_UCell, fill=cell_type)) + 
  geom_violin() +  geom_boxplot(width=0.1)

p

sodc <- AddModuleScore(sodc, features = list(mids), name=module)
sodcmeta <- sodc@meta.data
FeaturePlot(sodc, features = paste(module, "1", sep=""), label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  labs(title = paste(toupper(module), "enrichment"))

p <- ggplot(sodcmeta, aes(x=cell_type, y=m11, fill=cell_type)) + 
  geom_violin() +  geom_boxplot(width=0.1)

p

ggsave(paste("../plots/final/", module, "_sc.png",sep=""), width = 7, height = 6, scale = 0.8, dpi = 150)

#######

a <- DotPlot(sonorm, features = mids)$data
#a <- DotPlot(sonorm, features = myeloid)$data
af <- a[a$pct.exp > 20,]
#m1 20 2
#m2 10 1
#m3 3 1
#m4 2 1
#m5 3 1
af <- af[which(af$avg.exp.scaled > 2),]
selected_features <- unique(af$features.plot)
selected_features <- as.character(selected_features)
labels <- mlist[which(mlist$hsapiens_homolog_ensembl_gene %in% selected_features),]
labels <- labels[-which(duplicated(labels$hsapiens_homolog_ensembl_gene)),]
labels <- labels[match(selected_features, labels$hsapiens_homolog_ensembl_gene),]
selected_features <- factor(selected_features, levels=selected_features)

DotPlot(sonorm, features = selected_features, dot.min = 0) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  scale_x_discrete(labels = labels$hsymbol) + 
  RotatedAxis()

#m1 18
#m2 12
#m3 11
#m4 10
#m5 12
ggsave(paste("../plots/final/", module, "_sc_dot.png",sep=""), width = 12, height = 4, scale = 1, dpi = 100)

######

DimPlot(sonorm, label = TRUE, repel = TRUE, alpha=1) + NoLegend() + labs(title = "Heart cell atlas of normal tissue")
ggsave("../plots/final/sc_huebner.png", width = 7, height = 6, scale = 0.8, dpi = 150)




########## Pressure heart

pso <- readRDS("../metadata/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells.rds")
psometa <- pso@meta.data
psosub <- subset(pso, subset = gender == "Female")
psosub <- subset(psosub, subset = region == "LV")
psoas <- subset(psosub, subset = condition == "AS")
psonorm <- subset(psosub, subset = dataset == "hca")
Idents(psonorm) <- "cell_type"
Idents(psoas) <- "cell_type"
ftrs <- pso@assays$RNA@meta.features

modules <- c("m1", "m2", "m3", "m4", "m5")
markers <- list(m1=c(), m2=c(), m3=c(), m4=c(), m5=c())

for (i in 1:length(modules)) {
  mlist <- merge(eval(parse(text = modules[i])), t2gh, by="ens_gene")
  mlist <- mlist[-which(mlist$hsapiens_homolog_ensembl_gene==""), ]
  if (length(which(duplicated(mlist$hsapiens_homolog_ensembl_gene)))>0) { 
    mids <- mlist$hsymbol[-which(duplicated(mlist$hsapiens_homolog_ensembl_gene))]
  } else {
    mids <- mlist$hsymbol
  }
  
  markers[[i]] <- mids
}


factor(psometa$dataset)

psonorm <- AddModuleScore_UCell(psonorm, features = markers)
signature.names <- paste0(names(markers), "_UCell")
VlnPlot(psonorm, features = signature.names, group.by = "cell_type", pt.size=0)

psoas <- AddModuleScore_UCell(psoas, features = markers)
signature.names <- paste0(names(markers), "_UCell")
VlnPlot(psoas, features = signature.names, group.by = "cell_type", pt.size=0)
