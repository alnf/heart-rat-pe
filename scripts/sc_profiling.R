library(Seurat)
library(ggplot2)
library(RColorBrewer)

#so <- readRDS("../metadata/huebner_heart_cell_atlas.rds")
#remove(so)
#gc()
#someta <- so@meta.data
#head(Features(so))
#ftrs <- so@assays$RNA@meta.features
#sonorm <- subset(so, subset = disease == "normal")
#sonorm <- subset(sonorm, subset = Region_x == "LV")

#Idents(sonorm) <- "cell_type"


module = "m5"
mlist <- merge(eval(parse(text = module)), t2gh, by="ens_gene")
mlist <- mlist[-which(mlist$hsapiens_homolog_ensembl_gene==""), ]
mids <- mlist$hsapiens_homolog_ensembl_gene[-which(duplicated(mlist$hsapiens_homolog_ensembl_gene))]

sonorm <- AddModuleScore(sonorm, features = list(mids), name=module)
FeaturePlot(sonorm, features = paste(module, "1", sep=""), label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  labs(title = paste(toupper(module), "enrichment"))
ggsave(paste("../plots/final/", module, "_sc.png",sep=""), width = 7, height = 6, scale = 0.8, dpi = 150)

#######

a <- DotPlot(sonorm, features = mids)$data
af <- a[a$pct.exp > 3,]
#m1 20 2
#m2 10 1
#m3 3 1
#m4 2 1
#m5 3 1
af <- af[which(af$avg.exp.scaled > 1),]
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
psonorm <- subset(psonorm, subset = Region_x == "LV")
psometa <- psonorm@meta.data

psometa$disease

#sonorm <- subset(so, subset = disease == "normal")
