library(pheatmap)
library(CEMiTool)
library(CeTF)
library(dplyr)

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

metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$Pheno == "np"), ]

t2gh <- biomaRt::getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), mart = ensembl)
t2gh <- dplyr::rename(t2gh, ens_gene = ensembl_gene_id, hsymbol = hsapiens_homolog_associated_gene_name)

#LV: PEpp vs SDpp
lfc = 0.58

pp <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
pp <- pp[which(pp$padj<0.05),]
pp <- pp[which(abs(pp$log2FoldChange)>lfc),]
pp <- merge(pp, t2gh, by="ens_gene")
pp <- pp[-which(pp$hsymbol==""),]

#LV: PEd21 vs SDd21
lfc = 0.58

d21 <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21.tsv", sep="\t", header = T, stringsAsFactors = F)
d21 <- d21[which(d21$padj<0.05),]
d21 <- d21[which(abs(d21$log2FoldChange)>lfc),]
d21 <- merge(d21, t2gh, by="ens_gene")
d21 <- d21[-which(d21$hsymbol==""),]

#Heatmap 

lcounts <- as.data.frame(lcounts) 
m <- lcounts[all_genes_final_strict, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- as.matrix(m[,2:16])

color = hue_pal()(length(levels(factor(metadata.left$Pheno))))
names(color) <- levels(factor(metadata.left$Pheno))
annoCol <- list(Pheno = color)

pheatmap(m, annotation = metadata.left[,c("Pheno"), drop=F], annotation_colors = annoCol,
         filename = "../plots/heatmap_545_genes.png", width = 10, height = 7, scale = "row",
         show_rownames = F, show_colnames = F, treeheight_row = 0) 

ha = HeatmapAnnotation(Pheno = metadata.left[,c("Pheno")])
Heatmap(t(scale(t(m))), col=viridis(100, option="A"), show_row_names = F, show_row_dend = F, show_column_names = F,
        top_annotation = ha)


## ssGSEA


## MA plots
pp <- pp[which(abs(pp$log2FoldChange)>1),]

degs <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
degs <- merge(degs, pp[c("ens_gene", "hsymbol")], by="ens_gene", all.x=T)

degs <- degs[which(!is.na(degs$log2FoldChange)),]

degs$degs <- NA
degs$degs[which(degs$log2FoldChange>0)] <- "up"
degs$degs[which(degs$log2FoldChange<0)] <- "down"
degs$degs[which(is.na(degs$hsymbol))] <- NA
degs$degs <- factor(degs$degs, levels = c("up", "down"))

MAplot(degs, 7.5, 1, 0.5)


d21 <- d21[which(abs(d21$log2FoldChange)>1),]

degs <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21.tsv", sep="\t", header = T, stringsAsFactors = F)
degs <- merge(degs, d21[c("ens_gene", "hsymbol")], by="ens_gene", all.x=T)

degs <- degs[which(!is.na(degs$log2FoldChange)),]
degs <- degs[which((degs$log2FoldChange)< 19),]

degs$degs <- NA
degs$degs[which(degs$log2FoldChange>0)] <- "up"
degs$degs[which(degs$log2FoldChange<0)] <- "down"
degs$degs[which(is.na(degs$hsymbol))] <- NA
degs$degs <- factor(degs$degs, levels = c("up", "down"))

MAplot(degs, 10, 1, 0.5)

# Line plots

metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left$Pheno <- factor(metadata.left$Pheno, levels = c("SDd21", "SDpp", "np", "PEd21", "PEpp"))

pp <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp_short.tsv", sep="\t", header = T, stringsAsFactors = F)
d21 <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21_short.tsv", sep="\t", header = T, stringsAsFactors = F)
#all_degs <- unique(c(pp$ens_gene, d21$ens_gene))
all_degs <- unique(c(pp$ens_gene))
all_degs <- d21[-which(d21$ens_gene %in% pp$ens_gene),]$ens_gene


clusters <- degPatterns(lcounts[all_degs, which(colnames(lcounts) %in% metadata.left$SampleNumber)], metadata = metadata.left, minc = 10,
                        time = "Pheno", plot = T, groupDifference = 1.5)

p <- remove_geom(clusters$plot, "GeomSmooth")
p
p + facet_wrap(~cluster, ncol=1)



lines <- clusters$normalized[, c("genes", "Region", "Pheno", "value", "cluster")]
lines$cluster[which(lines$cluster==4)] <- 1
lines$cluster[which(lines$cluster==3)] <- 2

lines$cluster[which(lines$cluster==2)] <- 3
lines$cluster[which(lines$cluster==5)] <- 3


lines$phenotype <- "SD"
lines$phenotype[which(grepl("PE",lines$Pheno))] <- "PE"


ggplot(data=lines, aes(x=Pheno, y=value, group=genes, color=phenotype)) +
  geom_line(alpha=0.2, color="gray")+
  geom_point(position = position_jitterdodge(dodge.width = 0.1), alpha=0.5)+
  geom_boxplot(aes(group=Pheno, color=phenotype), alpha=0.5) +
  facet_wrap(~cluster)



ifng <- read.table("../degs/pathways/IFN_GAMMA_pe.tsv", sep="\t", header = F, stringsAsFactors = F)
colnames(ifng)[1] <- "hsymbol"
ifng <- merge(ifng, t2gh, by="hsymbol")
lines$ifng <- lines$genes %in% ifng$ens_gene 

ggplot(data=lines, aes(x=Pheno, y=value, group=genes, color=ifng)) +
  geom_line(alpha=0.2)+
  geom_point(position = position_jitterdodge(dodge.width = 0.5), alpha=0.5)+
  facet_wrap(~cluster)


lines <- clusters$normalized[, c("genes", "Region", "Pheno", "value", "cluster")]
inf_res <- read.table("../degs/pathways/INFLAMMATORY_RESPONSE_pe.tsv", sep="\t", header = F, stringsAsFactors = F)
colnames(inf_res)[1] <- "hsymbol"
inf_res <- merge(inf_res, t2gh, by="hsymbol")
lines$inf_res <- lines$genes %in% inf_res$ens_gene 

ggplot(data=lines, aes(x=Pheno, y=value, group=genes, color=inf_res)) +
  geom_line(alpha=0.2)+
  geom_point(position = position_jitterdodge(dodge.width = 0.5), alpha=0.5)+
  facet_wrap(~cluster)


# CEMiTool
metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$PhenoNames == "WTnp"), ]
m <- lcounts[, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m <- as.data.frame(m)
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
m <- m[m$hsymbol!="",]
rownames(m) <- m$hsymbol
m <- as.matrix(m[,2:(ncol(m)-1)])


anno_data <- metadata.left[,c("SampleNumber", "PhenoNames")]
anno_data$SampleNumber <- as.character(anno_data$SampleNumber)
colnames(anno_data) <- c("SampleName", "Class")
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
cem <- cemitool(as.data.frame(m), anno_data, filter = T)


cem
diagnostic_report(cem, force = T)
generate_report(cem, force = T)

gmt_in <- read_gmt("../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt")
cem <- mod_ora(cem, gmt_in)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]
plots[3]
plots[5]

gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
cem <- mod_ora(cem, gmt_in)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]
plots[2]
plots[3]
plots[4]
plots[5]

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
int_df$gene1symbol <- toupper(int_df$gene1symbol)
int_df$gene2symbol <- toupper(int_df$gene2symbol)
head(int_df)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem, n=20) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[3]
generate_report(cem, force = T)

idx <-  which(cem@module[which(cem@module$modules=="M3"),]$genes %in% pp$hsymbol)
cem@module[which(cem@module$modules=="M3"),]$genes[idx]


# AGT

m <- exprs.h[, which(colnames(exprs.h) %in% metadata.h$SampleNumber)]
m <- as.data.frame(m)
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
m <- m[m$hsymbol!="",]
rownames(m) <- m$hsymbol
m <- as.matrix(m[,2:(ncol(m)-1)])

AGT <- m[which(rownames(m)=="AGT"),]
AGT
