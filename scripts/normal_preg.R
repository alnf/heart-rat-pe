library(ssGSEA2)
library(cmapR)
library(pheatmap)

remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}

t2gh <- biomaRt::getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), mart = ensembl)
t2gh <- dplyr::rename(t2gh, ens_gene = ensembl_gene_id, hsymbol = hsapiens_homolog_associated_gene_name)

lfc = 0.58

#LV: SDpp vs np
ltEP <- read.table("../degs/WT/LV/degs_LV_SDpp_LV_np.tsv", sep="\t", header = T, stringsAsFactors = F)
ltEP <- ltEP[which(ltEP$padj<0.05),]
ltEP <- ltEP[which(abs(ltEP$log2FoldChange)>lfc),]
#ltEP <- read.table("../degs/WT/LV/degs_LV_SDpp_LV_np_short.tsv", sep="\t", header = T, stringsAsFactors = F)
ltEP <- merge(ltEP, t2gh, by="ens_gene")
ltEP <- ltEP[-which(ltEP$hsymbol==""),]

#LV: SDd21 vs np

RDP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_np.tsv", sep="\t", header = T, stringsAsFactors = F)
RDP <- RDP[which(RDP$padj<0.05),]
RDP <- RDP[which(abs(RDP$log2FoldChange)>lfc),]
#RDP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_np_short.tsv", sep="\t", header = T, stringsAsFactors = F)
RDP <- merge(RDP, t2gh, by="ens_gene")
RDP <- RDP[-which(RDP$hsymbol==""),]

#LV: SDd21 vs SDpp

RAP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
RAP <- RAP[which(RAP$padj<0.05),]
RAP <- RAP[which(abs(RAP$log2FoldChange)>lfc),]
#RAP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_SDpp_short.tsv", sep="\t", header = T, stringsAsFactors = F)
RAP <- merge(RAP, t2gh, by="ens_gene")
RAP <- RAP[-which(RAP$hsymbol==""),]

#Venn diagram

all_genes <- unique(c(RDP$hsymbol, RAP$hsymbol, ltEP$hsymbol))
fname = paste(c("RDP", "RAP", "ltEP"), collapse="_")
title = paste(fname, ", ngenes=", length(all_genes), sep="")
makeVenn(3, list(RDP$hsymbol, RAP$hsymbol, ltEP$hsymbol), c("RDP", "RAP", "ltEP"), title,
         paste("../plots/venn_regions_lfc1_", fname, ".png", sep=""),
         brewer.pal(3, "Pastel2"), dist = -0.35)

# Line plots

metadata.left <- metadata.h[which(grepl("SD",metadata.h$Pheno) | grepl("np",metadata.h$Pheno)), ]
metadata.left <- metadata.left[which(metadata.left$Region == "LV"), ]
metadata.left$Pheno <- as.character(metadata.left$Pheno)
metadata.left$Pheno <- factor(metadata.left$Pheno, levels = c("np", "SDd21", "SDpp"))

all_degs <- unique(c(RDP$ens_gene, RAP$ens_gene, ltEP$ens_gene))

clusters <- degPatterns(lcounts[all_degs, which(colnames(lcounts) %in% metadata.left$SampleNumber)], metadata = metadata.left, minc = 10,
                        time = "Pheno", plot = T, cutoff = 0.2)

p <- remove_geom(clusters$plot, "GeomSmooth")
p
p + facet_wrap(~cluster, ncol=1)

clus <- data.frame(ens_gene=clusters$df$genes, cluster=clusters$df$cluster)

# Heatmap

m <- lcounts[all_degs, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m <- as.data.frame(m)
m$ens_gene <- rownames(m)
m <- merge(m, clus, by="ens_gene")
m$RDP <- ""
m$RDP[which(m$ens_gene %in% RDP$ens_gene)] <- "True"
m$RAP <- ""
m$RAP[which(m$ens_gene %in% RAP$ens_gene)] <- "True"
m$ltEP <- ""
m$ltEP[which(m$ens_gene %in% ltEP$ens_gene)] <- "True"
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
cdesc <- data.frame(id=metadata.left$SampleNumber, type=as.character(metadata.left$Pheno))
rdesc <- data.frame(id=rownames(m), group=paste("group", m$cluster, sep="_"), RDP=m$RDP, RAP=m$RAP, ltEP=m$ltEP)
m <- as.matrix(m[,2:19])
(gct.data <- new("GCT", mat=m, cdesc=cdesc, rdesc=rdesc))

write_gct(gct.data, "../degs/norm_degs_058.gct")

res = run_ssGSEA2("../degs/norm_degs_058.gct_n18x546.gct",
                  output.prefix = "norm",
                  gene.set.databases = "../metadata/msigdb/h.all.v2023.2.Hs.symbols.gmt",
                  output.directory = "../degs/ssGSEA/",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 1000, 
                  min.overlap = 15, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "../degs/ssGSEA/run.log",
                  spare.cores = 20,
                  par=T)




# Heatmap
lcounts <- as.data.frame(lcounts) 
m <- lcounts[all_degs, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- as.matrix(m[,2:19])

color = hue_pal()(length(levels(factor(metadata.left$Pheno))))
names(color) <- levels(factor(metadata.left$Pheno))
annoCol <- list(Pheno = color)

pheatmap(m, annotation = metadata.left[,c("Pheno"), drop=F], annotation_colors = annoCol,
         filename = "../plots/heatmap_500_genes.png", width = 10, height = 7, scale = "row",
         main = "500 genes") 
