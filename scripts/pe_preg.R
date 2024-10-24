library(CEMiTool)
library(CeTF)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(viridis)
library(ggrepel)

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

clr <- read.table("../metadata/pe_colors.tsv", sep="\t", check.names = F, header = T, comment.char = "")
colors <- clr$color
names(colors) <- clr$variable
annoCol <- list(Group = colors)

metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$Pheno == "np"), ]

t2gh <- biomaRt::getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), mart = ensembl)
t2gh <- dplyr::rename(t2gh, ens_gene = ensembl_gene_id, hsymbol = hsapiens_homolog_associated_gene_name)

lcounts <- as.data.frame(lcounts)
agt <- t2gh[which(t2gh$hsymbol == "AGT"),]
lcounts <- lcounts[-which(rownames(lcounts)==agt$ens_gene),]
tpm_abd <- tpm_abd[-which(rownames(tpm_abd)==agt$ens_gene),]

#LV: PEpp vs SDpp
lfc = 1

pp <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
pp <- pp[-which(pp$ens_gene==agt$ens_gene),]
pp <- pp[which(pp$padj<0.05),]
pp <- pp[which(abs(pp$log2FoldChange)>lfc),]
pp <- merge(pp, t2gh, by="ens_gene")
pp <- pp[-which(pp$hsymbol==""),]
pp <- pp[-which(pp$symbol==""),]

#LV: PEd21 vs SDd21
lfc = 1

d21 <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21.tsv", sep="\t", header = T, stringsAsFactors = F)
d21 <- d21[-which(d21$ens_gene==agt$ens_gene),]
d21 <- d21[which(d21$padj<0.05),]
d21 <- d21[which(abs(d21$log2FoldChange)>lfc),]
d21 <- merge(d21, t2gh, by="ens_gene")
d21 <- d21[-which(d21$hsymbol==""),]
d21 <- d21[-which(d21$symbol==""),]

#Heatmap: all samples

m <- lcounts[all_genes_final_strict$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- as.matrix(m[,2:(ncol(m)-1)])
c <- cor(m)

cha = HeatmapAnnotation(Group = as.factor(metadata.left[,c("PhenoNames")]), col=annoCol, show_legend = F, show_annotation_name = F)
rha = rowAnnotation(Group = as.factor(metadata.left[,c("PhenoNames")]), col=annoCol, show_legend = F, show_annotation_name = F)

pals = hcl.pals("diverging")
pals
colors = hcl.colors(10, palette = "Purple-Green")
ht = Heatmap(c, show_row_names = T, show_row_dend = F, show_column_names = T,
        top_annotation = cha, right_annotation = rha, name = "corr", column_split = metadata.left[,c("PhenoNames")],
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 14),
        column_dend_height=unit(10, "mm"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14))
        )
png("../plots/pe_preg/heatmap_samples.png", width = 10.5, height = 9, units="in", res=80)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

#Heatmap: main genes

m <- lcounts[pp[which(abs(pp$log2FoldChange)>1),]$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

all_genes_final_strict <- all_genes_final_strict[-which(all_genes_final_strict$ens_gene==agt$ens_gene),,drop=F]
m <- lcounts[all_genes_final_strict$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

all_genes_WT_PE <- all_genes_WT_PE[-which(all_genes_WT_PE$ens_gene==agt$ens_gene),,drop=F]
m <- lcounts[all_genes_WT_PE$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- m[-which(m$hsymbol=="KLRC2"),]
m <- m[-which(m$hsymbol=="KLRC3"),]
m <- m[-which(m$hsymbol=="KLRC4"),]
m <- as.matrix(m[,2:(ncol(m)-1)])

m <- merge(m, t2g, by="ens_gene")
m <- m[!duplicated(m$symbol),]
m <- m[which(m$symbol!=""),]
rownames(m) <- m$symbol
m <- m[-grep("ribosomal RNA", m$description),]
m <- as.matrix(m[,2:(ncol(m)-4)])


m <- m[,match(metadata.left$SampleNumber, colnames(m))]
ha = HeatmapAnnotation(Group = metadata.left[,c("PhenoNames")], col=annoCol, show_legend = F, show_annotation_name = F)

d21 <- d21[which(d21$ens_gene %in% m$ens_gene),]
m <- m[!duplicated(m$ens_gene),]
d21 <- d21[order(d21$log2FoldChange, decreasing = T),]
m <- m[match(d21$symbol, rownames(m)),]


colors = hcl.colors(30, palette = "Purple-Green")
ht = Heatmap(t(scale(t(m))), show_row_names = T, show_row_dend = F, show_column_names = T, cluster_columns = F,
             top_annotation = ha, name = "expr", column_split = cspl)
ht = Heatmap(t(scale(t(m))), show_row_names = F, show_row_dend = F, show_column_names = T, cluster_columns = T, cluster_rows = F,
             top_annotation = ha, name = "expr", column_split = cspl)

png("../plots/pe_preg/heatmap_main_genes.png", width = 10.5, height = 19, units="in", res=80)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

#Heatmap: all 600 genes

m <- lcounts[all_genes_final_strict$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- m[m$hsymbol!="",]
m <- as.matrix(m[,2:(ncol(m)-1)])

mod_cols = c(M1="#820D3F",M2="#E64A00",M3="#3B3EDE",M4="#871C9A",M5="#14C7BA", Not.Correlated="gray")

genes <- cem@module[which(cem@module$modules=="M2"),]$genes
m2 <- rep(NA, nrow(m))
names(m2) <- rownames(m)
m2[which(names(m2) %in% genes)] <- mod_cols[["M2"]]
genes <- cem@module[which(cem@module$modules=="M3"),]$genes
m3 <- rep(NA, nrow(m))
names(m3) <- rownames(m)
m3[which(names(m3) %in% genes)] <- mod_cols[["M3"]]
genes <- cem@module[which(cem@module$modules=="M4"),]$genes
m4 <- rep(NA, nrow(m))
names(m4) <- rownames(m)
m4[which(names(m4) %in% genes)] <- mod_cols[["M4"]]
genes <- cem@module[which(cem@module$modules=="M5"),]$genes
m5 <- rep(NA, nrow(m))
names(m5) <- rownames(m)
m5[which(names(m5) %in% genes)] <- mod_cols[["M5"]]
genes <- cem@module[which(cem@module$modules=="M1"),]$genes
m1 <- rep(NA, nrow(m))
names(m1) <- rownames(m)
m1[which(names(m1) %in% genes)] <- mod_cols[["M1"]]
    
rspl <- data.frame(modules = cem@module$modules)
rspl$hsymbol <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M5", "M3", "M4", "NC"))
#rspl <- rspl[order(rspl$modules),]
#m <- m[match(rspl$hsymbol, rownames(m)),]

annoRow <- list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5)
anno_df <- data.frame(M1 = rownames(m), M2 = rownames(m), M3 = rownames(m), M4 = rownames(m), M5 = rownames(m))
rha = rowAnnotation(df=anno_df, col=annoRow, show_legend = F)

metadata.left$PhenoNames <- factor(metadata.left$PhenoNames, levels = c("WTpreg", "WTpost", "PEpreg", "PEpost"))
metadata.left <- metadata.left[order(metadata.left$PhenoNames),]
m <- m[,match(metadata.left$SampleNumber, colnames(m))]

cspl <- metadata.left$PhenoNames

cha = HeatmapAnnotation(Group = metadata.left[,c("PhenoNames")], col=annoCol, show_legend = F, show_annotation_name = F)

colors = hcl.colors(30, palette = "Purple-Green")
ht = Heatmap(t(scale(t(m))), show_row_names = F, show_row_dend = F, show_column_names = T, cluster_columns = F,
             top_annotation = cha, right_annotation = rha, name = "expr",
             column_split = cspl, row_split = rspl$modules, cluster_row_slices = F,
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 14),
             column_dend_height=unit(10, "mm"),
             heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14))
             )
png("../plots/pe_preg/heatmap_all_genes.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

## MA plots
imp <- pp[which(abs(pp$log2FoldChange)>1),]

degs <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
degs <- merge(degs, imp[c("ens_gene", "hsymbol")], by="ens_gene", all.x=T)
degs <- degs[which(!is.na(degs$log2FoldChange)),]

degs$degs <- NA
degs$degs[which(degs$log2FoldChange>0)] <- "up"
degs$degs[which(degs$log2FoldChange<0)] <- "down"
degs$degs[which(is.na(degs$hsymbol))] <- NA
degs$degs <- factor(degs$degs, levels = c("up", "down"))

MAplot(degs, 7.5, 1, 0.5, "PEpost vs WTpost, logFC>1, ngenes = 139")
ggsave("../plots/pe_preg/ma_pp.png", width = 10.5, height = 6, dpi = 80)

imp <- d21[which(abs(d21$log2FoldChange)>1),]

degs <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21.tsv", sep="\t", header = T, stringsAsFactors = F)
degs <- merge(degs, d21[c("ens_gene", "hsymbol")], by="ens_gene", all.x=T)

degs <- degs[which(!is.na(degs$log2FoldChange)),]
degs <- degs[which((degs$log2FoldChange)< 19),]

degs$degs <- NA
degs$degs[which(degs$log2FoldChange>0)] <- "up"
degs$degs[which(degs$log2FoldChange<0)] <- "down"
degs$degs[which(is.na(degs$hsymbol))] <- NA
degs$degs <- factor(degs$degs, levels = c("up", "down"))

MAplot(degs, 10, 1, 0.5, "PEpreg vs WTpreg, logFC>1, ngenes = 510")
ggsave("../plots/pe_preg/ma_d21.png", width = 10.5, height = 6, dpi = 80)


## ssGSEA


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
m <- lcounts[all_genes_final$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m <- lcounts[all_genes_final_strict$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m <- lcounts[, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

all_genes_WT_PE <- all_genes_WT_PE[-which(all_genes_WT_PE$ens_gene==agt$ens_gene),,drop=F]
m <- lcounts[all_genes_WT_PE$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

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
cem <- cemitool(as.data.frame(m), anno_data, filter = F, filter_pval = 0.1,
                force_beta = T, apply_vst = F, rank_method = "median",
                gsea_min_size = 10
                )
cem
cem@enrichment_plot

gmt_in <- read_gmt("../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c5.go.bp.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c2.cp.biocarta.v2023.2.Hs.symbols.gmt")


gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

cem <- mod_ora(cem, gmt_in)
mod_colors(cem) <- mod_cols
cem <- plot_ora(cem, pv_cut=0.05)
plots <- show_plot(cem, "ora")

ora_bp <- cem@ora[which(cem@ora$Module=="M2"),]

plots[1]
ggsave("../plots/pe_preg/m1_enrichplot.png", width = 6, height = 5, dpi = 100)
plots[2]
ggsave("../plots/pe_preg/m2_enrichplot_bp.png", width = 7, height = 5, dpi = 100)
plots[3]
ggsave("../plots/pe_preg/m3_enrichplot.png", width = 6, height = 5, dpi = 100)
plots[4]
ggsave("../plots/pe_preg/m4_enrichplot.png", width = 6, height = 5, dpi = 100)
plots[5]
ggsave("../plots/pe_preg/m5_enrichplot.png", width = 6, height = 5, dpi = 100)

epl <- cem@enrichment_plot$enrichment_plot
epl + theme(axis.text.y = element_text(colour = mod_cols[c(1,3,4,2,5)], size=14, face = "bold"),
            axis.text.x = element_text(colour = annoCol$Group[c(4,3,2,1)], size=14, face = "bold"),
            axis.title = element_text(size=13)) +
      labs(x = "Group")
ggsave("../plots/pe_preg/main_enrichplot.png", width = 9, height = 6, dpi = 100)

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
int_df$gene1symbol <- toupper(int_df$gene1symbol)
int_df$gene2symbol <- toupper(int_df$gene2symbol)
head(int_df)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem, n=10) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]$M1 + theme(title = element_blank())
ggsave("../plots/pe_preg/m1_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[2]$M2 + theme(title = element_blank())
ggsave("../plots/pe_preg/m2_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[3]$M3 + theme(title = element_blank())
ggsave("../plots/pe_preg/m3_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[4]$M4 + theme(title = element_blank())
ggsave("../plots/pe_preg/m4_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[5]$M5 + theme(title = element_blank())
ggsave("../plots/pe_preg/m5_interactionplot.png", width = 7, height = 5, dpi = 80)


#generate_report(cem, force = T, directory="../degs/CEMi/degs_hallmark", title = "degs_hallmark")
generate_report(cem, force = T, directory="../degs/CEMi/degs_bp", title = "degs_bp")

idx <-  which(cem@module[which(cem@module$modules=="M5"),]$genes %in% pp$hsymbol)
m2g <- cem@module[which(cem@module$modules=="M2"),]$genes
write.table(m2g, "../degs/module2.tsv", sep="\t", row.names = F)

# Module heatmap

genes <- cem@module[which(cem@module$modules=="M5"),]$genes
mat <- m[genes,]

ha = HeatmapAnnotation(Group = metadata.left[,c("PhenoNames")], col=annoCol)
colors = hcl.colors(30, palette = "Purple-Green")
ht = Heatmap(t(scale(t(mat))), show_row_names = T, col = colors, show_row_dend = F, show_column_names = T,
             top_annotation = ha, name = "expr",)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))

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


enrichPlot()
cem@enrichment_plot$enrichment_plot$theme$text$size = 13
cem@enrichment_plot

cem@ora$geneID[which(cem@ora$Module=="M2")]

# performing getEnrich analysis
pp <- pp[-which(pp$hsymbol=="KLRC2"),]
pp <- pp[-which(pp$hsymbol=="KLRC3"),]
pp <- pp[-which(pp$hsymbol=="KLRC4"),]
pp <- pp[-which(duplicated(pp$hsymbol)),]

pp_cond <- getEnrich(genes = pp$symbol, organismDB = org.Rn.eg.db, keyType = 'SYMBOL', 
                   ont = 'BP', fdrMethod = "BH", fdrThr = 0.05, minGSSize = 10, 
                   maxGSSize = 500)

d21_cond <- getEnrich(genes = d21$symbol, organismDB = org.Rn.eg.db, keyType = 'SYMBOL', 
                     ont = 'BP', fdrMethod = "BH", fdrThr = 0.05, minGSSize = 10, 
                     maxGSSize = 500)

pp_cond_res <- pp_cond$results
pp_cond_res <- pp_cond_res[pp_cond_res$p.adjust<=0.05,]


p <- enrichPlot(res = pp_cond$results, showCategory = 10,
                type = "circle") +
                theme(axis.text.x = element_text(size=12))
p

enrichPlot(res = d21_cond$results, showCategory = 10,
           type = "circle") +
           theme(axis.text.x = element_text(size=12))

d21_cond_res <- d21_cond$results
d21_cond_res <- d21_cond_res[d21_cond_res$p.adjust<=0.05,]


### Venn diagram: SD vs PE
region = "LV"
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/pe_preg/venn_", fname, ".png", sep="")
names = c("PEpreg_SDpreg", "PEpost_SDpost")
makeVenn(2, list(d21$hsymbol, pp$hsymbol), names, "DEGs in pregnancy and postpartum",
         fname, c("lightgreen", "pink"))


duplicated(intersect(d21$hsymbol, pp$hsymbol))


library(ggVennDiagram)

# List of items
x <- list(A = 1:5, B = 2:7)

# 2D Venn diagram
ggVennDiagram(list(d21$ens_gene, pp$ens_gene))

duplicated(pp$hsymbol)


# CEM for pp only

anno_data <- metadata.left[,c("SampleNumber", "PhenoNames")]
anno_data$SampleNumber <- as.character(anno_data$SampleNumber)
colnames(anno_data) <- c("SampleName", "Class")
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

ppp <- pp[which(pp$hsymbol %in% rownames(m)),]
mm <- m[ppp[which(ppp$log2FoldChange>0),]$hsymbol,]
kk <- m[ppp[which(ppp$log2FoldChange<0),]$hsymbol,]



cem <- cemitool(as.data.frame(kk), anno_data, filter = F, filter_pval = 0.1,
                force_beta = T, apply_vst = F, rank_method = "median",
                gsea_min_size = 10
)
cem
cem@enrichment_plot


gmt_in <- read_gmt("../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c5.go.bp.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c2.cp.biocarta.v2023.2.Hs.symbols.gmt")

gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

cem <- mod_ora(cem, gmt_in)
#mod_colors(cem) <- mod_cols
cem <- plot_ora(cem, pv_cut=0.5)
plots <- show_plot(cem, "ora")


#cem@enrichment_plot
plots[1]

ora_selected <- cem@ora





# CEM for final heatmaps

anno_data <- metadata.left[,c("SampleNumber", "PhenoNames")]
anno_data$SampleNumber <- as.character(anno_data$SampleNumber)
colnames(anno_data) <- c("SampleName", "Class")
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

m <- lcounts[genes_pe_lfc1$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m <- lcounts[genes_pe_reduced_lfc1$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

cem <- cemitool(as.data.frame(m), anno_data, filter = F,
                force_beta = T, apply_vst = F, rank_method = "median",
                gsea_min_size = 10
)
cem
cem@enrichment_plot

m2genes <- cem@module$genes[which(cem@module$modules=="M2")]
m2genes <- m2genes[which(m2genes %in% d21$ens_gene)]
lfcsign <- ifelse(d21$log2FoldChange[d21$ens_gene %in% m2genes] > 0, "P1", "P2")
module <- data.frame(genes=m2genes, modules=lfcsign)
cem@module <- module

gmt_in <- read_gmt("../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c5.go.bp.v2023.2.Hs.symbols.gmt")
gmt_in <- read_gmt("../databases/msigdb/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")
gmt_in <- merge(gmt_in, t2gh, by.x="gene", by.y="hsymbol")

gmt_merged <- gmt_in
gmt_in <- gmt_in[,c(2:3)]
colnames(gmt_in)[2] <- "gene" 

cem <- mod_ora(cem, gmt_in)
#mod_colors(cem) <- mod_cols
cem <- plot_ora(cem, pv_cut=0.05)
plots <- show_plot(cem, "ora")
plots[1]
plots[2]
plots[3]
plots[4]
plots[5]

figure <- ggarrange(plots[1]$M1$pl, plots[2]$M2$pl, plots[3]$M3$pl, plots[4]$M4$pl, plots[5]$M5$pl,
                    ncol = 3, nrow = 2)
figure <- ggarrange(plots[1]$M1$pl, plots[2]$M2$pl, plots[3]$M3$pl,
                    ncol = 3, nrow = 1)
figure

ora <- cem@ora
cem@module$genes[which(cem@module$modules=="M4")]

infgamma <- ora[which(ora$ID=="HALLMARK_INTERFERON_GAMMA_RESPONSE"),]$geneID
infgamma <- str_split(infgamma, "/")
infgamma <- unlist(infgamma)
infgamma <- unique(infgamma)

infalpha <- ora[which(ora$ID=="HALLMARK_INTERFERON_ALPHA_RESPONSE"),]$geneID
infalpha <- str_split(infalpha, "/")
infalpha <- unlist(infalpha)
infalpha <- unique(infalpha)
infs <- c(infgamma, infalpha)
infs <- unique(infs)

e2f <- ora[which(ora$ID=="HALLMARK_E2F_TARGETS"),]$geneID
e2f <- str_split(e2f, "/")
e2f <- unlist(e2f)
e2f <- unique(e2f)

mis <- ora[which(ora$ID=="HALLMARK_MITOTIC_SPINDLE"),]$geneID
mis <- str_split(mis, "/")
mis <- unlist(mis)
mis <- unique(mis)

g2m <- ora[which(ora$ID=="HALLMARK_G2M_CHECKPOINT"),]$geneID
g2m <- str_split(g2m, "/")
g2m <- unlist(g2m)
g2m <- unique(g2m)

gmt_merged[gmt_merged$ens_gene %in% infs,]
gmt_merged[gmt_merged$ens_gene %in% cem@module$genes[which(cem@module$modules=="M4")],]


ora[which(ora$Module=="M4"),]$geneID


ora_all <- cem@ora
ora_back <- ora_all
ora_all <- rbind(ora_all, cem@ora)
ora_all <- ora_all[which(ora_all$p.adjust<0.05),]
ora_all <- ora_all[which(!is.na(ora_all$qvalue)),]

ora_all$pathway <- sapply(strsplit(ora_all$ID,"_"), `[`, 1)
ora_all$ID <- sub(".*?_", "", ora_all$ID)
ora_all$ID <- gsub("_", " ", ora_all$ID)
ora_all$geneRatio <- sapply(ora_all$GeneRatio, function(x) eval(parse(text=x)))

ora_all <- ora_all[which(ora_all$Count>5),]

ora_all %>% group_by(Module) %>% slice_max(order_by = geneRatio, n = 10) -> ora_all
ora_all$pathway <- as.factor(ora_all$pathway)
ora_all$ID <- factor(ora_all$ID, levels=ora_all$ID)

ordered_ids <- forcats::fct_reorder(ora_all$ID, ora_all$geneRatio)
ora_all <- ora_all[order(match(ora_all$ID, ordered_ids)),]

ora_all %>% group_by(Module)

a <- ifelse(ora_all$pathway == "HALLMARK", "black", "#217364")
names(a) <- ora_all$pathway

library(ggh4x)
strip <- strip_themed(text_x = elem_list_text(color = mcols))

mcols <- c(M1="#820D3F", M2="#E64A00", M3="#3B3EDE", M4="#871C9A", M5="#14C7BA")
png("../plots/enrich_genes_pe_lfc1.png",  width = 17, height = 6, units="in", res=80)

png("../plots/enrich_genes_pe_reduced_lfc1.png",  width = 20, height = 5, units="in", res=80)
ggplot(ora_all,
       aes(geneRatio, forcats::fct_reorder(ID, geneRatio))) + 
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
       scale_size_continuous(range=c(1, 7)) +
       geom_segment(aes(xend=0, yend = ID)) +
       theme_minimal() +
       theme(axis.text.y = element_text(colour = a)) +
       xlab("Gene Ratio") +
       ylab(NULL) + 
       facet_wrap2(~Module, scales="free", strip = strip) +
       theme(strip.text.x = element_text(colour = mcols, size=14, face="bold"))

pcols <- c("HALLMARK"="black", "GOBP"="#217364")
lgd = Legend(labels = names(pcols), title = "Pathway", labels_gp = gpar(fontsize = 8), nrow = 1, legend_gp = gpar(fill = pcols))
draw(lgd, x = unit(0.3, "in"), y = unit(0.5, "in"), just = c("left", "bottom"))
dev.off()


ggplot(ora_all,
       aes(geneRatio, forcats::fct_reorder(ID, geneRatio))) + 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  geom_segment(aes(xend=0, yend = ID)) +
  theme_minimal() +
  theme(axis.text.y = element_text(colour = a)) +
  xlab("Gene Ratio") +
  ylab(NULL) + 
  facet_wrap(~Module, scales="free") +
  theme(strip.text.x = element_text(colour = mcols, size=14, face="bold"))



# Upset plot
ora_back <- ora_all
ora_all$geneID <- sapply(strsplit(ora_all$geneID,"/"), `[`, )
names(ora_all$geneID) <- ora_all$ID

modname = "M2"
ups <- make_comb_mat(ora_all[which(ora_all$Module==modname),]$geneID, mode = "intersect")
ups <- ups[comb_size(ups) >= 0 & comb_degree(ups) == 2]

us <- UpSet(ups, pt_size = unit(4, "mm"), lwd = 2,
      comb_col = brewer.pal(comb_degree(ups)[1], "Dark2")[comb_degree(ups)],
      top_annotation = upset_top_annotation(ups, add_numbers = TRUE),
      right_annotation = upset_right_annotation(ups, add_numbers = TRUE),
      row_names_gp = grid::gpar(fontsize = 9),
      comb_order = order(comb_size(ups)))
#png(paste("../plots/final/final_upset_genes_pe_lfc1_",modname,".png", sep=""), width = 15, height = 5, units="in", res=100)
draw(us, padding = unit(c(2, 29, 2, 2), "mm"))
#dev.off()

library(dbscan)
library(ggrepel)
library("FactoMineR")
library("factoextra")

oras <- ora_all[which(ora_all$Module==modname),]
ltm <- list_to_matrix(oras$geneID)
ltm <- as.data.frame(ltm)
ltm[] <- lapply(ltm, as.character)
cats = apply(ltm, 2, function(x) nlevels(as.factor(x)))
cats

mca <- MCA(t(ltm), ncp = 2, graph = F)
fviz_mca_ind(mca, 
             repel = TRUE, # Avoid text overlapping (slow if many point)
             ggtheme = theme_minimal())


mca_df = data.frame(mca$var$coord, Variable = rep(names(cats), cats))
autoplot(mca_df, data = oras, colour = "Cluster",
         label.show.legend = F) +
  geom_text_repel(aes(label = ID, colour = Cluster), size = 3, box.padding = 1.5)




ltm <- list_to_matrix(oras$geneID)
pca <- prcomp(t(ltm))
rownames(oras) <- oras$ID
cl <- hdbscan(t(ltm), minPts = 2, gen_simplified_tree = F)
oras$Cluster <- as.factor(cl$cluster)
autoplot(pca, data = oras, colour = "Cluster",
         label.show.legend = F) +
  geom_text_repel(aes(label = ID, colour = Cluster), size = 3, max.overlaps = 25, force = 2)



