library(ssGSEA2)
library(cmapR)
library(pheatmap)
library(ggpubr)
library(ComplexHeatmap)
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

t2gh <- biomaRt::getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), mart = ensembl)
t2gh <- dplyr::rename(t2gh, ens_gene = ensembl_gene_id, hsymbol = hsapiens_homolog_associated_gene_name)

lfc = 0.58

#LV: SDpp vs np
ltEP <- read.table("../degs/LV/degs_LV_SDpp_LV_np.tsv", sep="\t", header = T, stringsAsFactors = F)
ltEP <- ltEP[which(ltEP$padj<0.05),]
ltEP <- ltEP[which(abs(ltEP$log2FoldChange)>lfc),]
#ltEP <- read.table("../degs/WT/LV/degs_LV_SDpp_LV_np_short.tsv", sep="\t", header = T, stringsAsFactors = F)
ltEP <- merge(ltEP, t2gh, by="ens_gene")
ltEP <- ltEP[-which(ltEP$hsymbol==""),]
#ltEP <- ltEP[-which(ltEP$ens_gene==agt$ens_gene),]

#LV: SDd21 vs np

RDP <- read.table("../degs/LV/degs_LV_SDd21_LV_np.tsv", sep="\t", header = T, stringsAsFactors = F)
RDP <- RDP[which(RDP$padj<0.05),]
RDP <- RDP[which(abs(RDP$log2FoldChange)>lfc),]
#RDP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_np_short.tsv", sep="\t", header = T, stringsAsFactors = F)
RDP <- merge(RDP, t2gh, by="ens_gene")
RDP <- RDP[-which(RDP$hsymbol==""),]
#RDP <- RDP[-which(RDP$ens_gene==agt$ens_gene),]

#LV: SDd21 vs SDpp

RAP <- read.table("../degs/LV/degs_LV_SDd21_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
RAP <- RAP[which(RAP$padj<0.05),]
RAP <- RAP[which(abs(RAP$log2FoldChange)>lfc),]
#RAP <- read.table("../degs/WT/LV/degs_LV_SDd21_LV_SDpp_short.tsv", sep="\t", header = T, stringsAsFactors = F)
RAP <- merge(RAP, t2gh, by="ens_gene")
RAP <- RAP[-which(RAP$hsymbol==""),]
#RAP <- RAP[-which(RAP$ens_gene==agt$ens_gene),]


RAP <- RAP[order(abs(RAP$log2FoldChange), decreasing = T),]
RAP$logMean <- log(RAP$baseMean)
RAP$rank <- abs(RAP$logMean*RAP$log2FoldChange)



#Venn diagram

all_genes <- unique(c(RDP$hsymbol, RAP$hsymbol, ltEP$hsymbol))
fname = paste(c("RDP", "RAP", "ltEP"), collapse="_")
title = paste(fname, ", ngenes=", length(all_genes), sep="")
makeVenn(3, list(RDP$hsymbol, RAP$hsymbol, ltEP$hsymbol), c("RDP", "RAP", "ltEP"), title,
         paste("../plots/normal_preg/venn_regions_lfc0.58_", fname, ".png", sep=""),
         brewer.pal(3, "Pastel2"), dist = -0.35)

# Line plots

metadata.left <- metadata.h[which(grepl("SD",metadata.h$Pheno) | grepl("np",metadata.h$Pheno)), ]
metadata.left <- metadata.left[which(metadata.left$Region == "LV"), ]
metadata.left$Pheno <- as.character(metadata.left$Pheno)
metadata.left$Pheno[which(metadata.left$Pheno=="SDd21")] <- "preg"
metadata.left$Pheno[which(metadata.left$Pheno=="SDpp")] <- "pp"

metadata.left$Pheno <- factor(metadata.left$Pheno, levels = c("np", "preg", "pp"))



all_degs <- unique(c(RDP$ens_gene, RAP$ens_gene, ltEP$ens_gene))

clusters <- degPatterns(lcounts[all_degs, which(colnames(lcounts) %in% metadata.left$SampleNumber)], metadata = metadata.left, minc = 10,
                        time = "Pheno", plot = T, cutoff = 0.2)

p <- remove_geom(clusters$plot, "GeomSmooth")
p
p + facet_wrap(~cluster, ncol=1)

clus <- data.frame(ens_gene=clusters$df$genes, cluster=clusters$df$cluster)
write.table(clus[which(clus$cluster==1),"ens_gene"], "../degs/normal_preg_cluster1.tsv", sep="\t", row.names = F, col.names = T)
write.table(clus[which(clus$cluster==2),"ens_gene"], "../degs/normal_preg_cluster2.tsv", sep="\t", row.names = F, col.names = T)
write.table(clus[which(clus$cluster==3),"ens_gene"], "../degs/normal_preg_cluster3.tsv", sep="\t", row.names = F, col.names = T)
write.table(clus[which(clus$cluster==4),"ens_gene"], "../degs/normal_preg_cluster4.tsv", sep="\t", row.names = F, col.names = T)

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
logFC <- ltEP$ens_gene[which(abs(ltEP$log2FoldChange) > 1)]
logFC <- c(logFC, RDP$ens_gene[which(abs(RDP$log2FoldChange) > 1)])
logFC <- c(logFC, RAP$ens_gene[which(abs(RAP$log2FoldChange) > 1.4)])
m$logFC <- ""
m$logFC[which(m$ens_gene %in% logFC)] <- "True"

m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
cdesc <- data.frame(id=metadata.left$SampleNumber, type=as.character(metadata.left$Pheno))
rdesc <- data.frame(id=rownames(m), group=paste("group", m$cluster, sep="_"), RDP=m$RDP, RAP=m$RAP, ltEP=m$ltEP, logFC=m$logFC)
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
         filename = "../plots/heatmap_545_genes.png", width = 10, height = 7, scale = "row",
         show_rownames = F, show_colnames = F, treeheight_row = 0) 

Heatmap(t(scale(t(m))), col=viridis(100, option="A"), show_row_names = F, show_row_dend = F, show_column_names = F)

# Line plot with pathways
lines <- clusters$normalized[, c("genes", "Region", "Pheno", "value", "cluster")]


ifna <- read.table("../degs/pathways/IFN_ALPHA_norma.tsv", sep="\t", header = F, stringsAsFactors = F)
colnames(ifna)[1] <- "hsymbol"
ifna <- merge(ifna, t2gh, by="hsymbol")
lines$ifna <- lines$genes %in% ifna$ens_gene 

spindle <- read.table("../degs/pathways/MITOTIC_SPINDLE_norma.tsv", sep="\t", header = F, stringsAsFactors = F)
colnames(spindle)[1] <- "hsymbol"
spindle <- merge(spindle, t2gh, by="hsymbol")
lines$spindle <- lines$genes %in% spindle$ens_gene 

  
ggplot(data=lines, aes(x=Pheno, y=value, group=genes, color=ifna)) +
  geom_line(alpha=0.2)+
  geom_point(position = position_jitterdodge(dodge.width = 0.9), alpha=0.5)+
  facet_wrap(~cluster)

ggplot(data=lines, aes(x=Pheno, y=value, group=genes, color=spindle)) +
  geom_line(alpha=0.2)+
  geom_point(position = position_jitterdodge(dodge.width = 0.9), alpha=0.5)+
  facet_wrap(~cluster)

table(lines$cluster)


stat <- read.table("../metadata/groups_line_plot.tsv", sep="\t", header = T, stringsAsFactors = F)
stat <- stat %>%
  mutate(y.position = c(1.5, 1.7, 1.9))

stat$p <- as.factor(stat$p)
stat$col <- stat$p

nclus = 4
ggplot(data=lines[which(lines$cluster==nclus),], aes(x=Pheno, y=value, group=genes, fill = "lightgray")) +
  geom_line(alpha=0.2, color = "gray")+
  geom_point(position = position_jitterdodge(dodge.width = 0.2), alpha=0.5) +
  stat_pvalue_manual(stat, label = "p", color = "col") +
  theme(legend.position="none") +
  xlab("Animal group") + ylab("z-score") +
  labs(title = paste("Expressional Pattern", nclus, sep=" "),
       subtitle = paste("ngenes:", length(levels(factor(lines[which(lines$cluster==nclus),]$genes))), sep=" "))

ggplot(data=lines[which(lines$cluster==nclus),], aes(x=Pheno, y=value, group=genes, fill = "lightgray")) +
  geom_line(alpha=0.2, color = "gray")+
  geom_point(position = position_jitterdodge(dodge.width = 0.2), alpha=0.5) +
  geom_bracket(
    xmin = "np", xmax = "d21", y.position = 1.7,
    label = "RDP"
  ) +
  theme(legend.position="none") +
  xlab("Animal group") + ylab("z-score") +
  labs(title = paste("Expressional Pattern", nclus, sep=" "),
       subtitle = paste("ngenes:", length(levels(factor(lines[which(lines$cluster==nclus),]$genes))), sep=" "))

  
palette.colors(palette = "Okabe-Ito")
palette.pals()
show_col(palette.colors(palette = "Okabe-Ito"))
show_col(palette.colors(palette = "Polychrome 36"))
show_col(palette.colors(palette = "Paired"))
show_col(palette.colors(palette = "R4"))


library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)

library(scales)
show_col(viridis(n = 8, option = "H"))
display.brewer.all(colorblindFriendly = TRUE)


library(colorBlindness)
cvdPlot(replacePlotColor(displayColors(palette.colors(palette = "Okabe-Ito"))))
cvdPlot(replacePlotColor(displayColors(viridis(n = 8, option = "H"))))

cvdPlot(replacePlotColor(show_col(palette.colors(palette = "Paired"))))

library(colorblindr)
p <- displayColors(palette.colors(palette = "Paired"))
p <- displayColors(palette.colors(palette = "Okabe-Ito"))
p <- displayColors(palette.colors(palette = "Set 1"))
cvd_grid(p)


# MA plot
degs <- read.table("../degs/LV/degs_LV_SDd21_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)
degs <- merge(degs, RAP[c("ens_gene", "hsymbol")], by="ens_gene", all.x=T)

degs <- degs[which(degs$ens_gene != "ENSRNOG00000032348"),]
degs <- degs[which(!is.na(degs$log2FoldChange)),]
degs <- degs[which((degs$log2FoldChange)> -19),]

degs$degs <- NA
degs$degs[which(degs$log2FoldChange>0)] <- "up"
degs$degs[which(degs$log2FoldChange<0)] <- "down"
degs$degs[which(is.na(degs$hsymbol))] <- NA
degs$degs <- factor(degs$degs, levels = c("up", "down"))

MAplot(degs, 10, 1, 0.5)


# CEMi

all_genes <- unique(c(RDP$ens_gene, RAP$ens_gene, ltEP$ens_gene))
m <- lcounts[all_genes, which(colnames(lcounts) %in% metadata.left$SampleNumber)]

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

#ont_type = "default"
#ont_type = "bp"
#ont_type = "hallmark"
#ont_type = "kegg"
ont_type = "biocarta"

plots[1]
ggsave(paste("../plots/normal_preg/m1_enrichplot_", ont_type, ".png", sep=""), width = 6, height = 5, dpi = 100)
plots[2]
ggsave(paste("../plots/normal_preg/m2_enrichplot_", ont_type, ".png", sep=""), width = 7, height = 5, dpi = 100)
plots[3]
ggsave(paste("../plots/normal_preg/m3_enrichplot_", ont_type, ".png", sep=""), width = 6, height = 5, dpi = 100)
plots[4]
ggsave(paste("../plots/normal_preg/m4_enrichplot_", ont_type, ".png", sep=""), width = 6, height = 5, dpi = 100)
plots[5]
ggsave(paste("../plots/normal_preg/m5_enrichplot_", ont_type, ".png", sep=""), width = 6, height = 5, dpi = 100)


annoColnp <- annoCol
annoColnp$Group <- c("WTnp"="#719360", "WTpost"="#CCBB44", "WTpreg"="#4477AA")

epl <- cem@enrichment_plot$enrichment_plot
epl + theme(axis.text.y = element_text(colour = mod_cols[c(1,3,4,2,5)], size=14, face = "bold"),
            axis.text.x = element_text(colour = annoColnp$Group[c(1,2,3)], size=14, face = "bold"),
            axis.title = element_text(size=13)) +
  labs(x = "Group")
ggsave("../plots/normal_preg/main_enrichplot.png", width = 7, height = 5, dpi = 100)

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
int_df$gene1symbol <- toupper(int_df$gene1symbol)
int_df$gene2symbol <- toupper(int_df$gene2symbol)
head(int_df)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem, n=10) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]$M1 + theme(title = element_blank())
ggsave("../plots/normal_preg/m1_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[2]$M2 + theme(title = element_blank())
ggsave("../plots/normal_preg/m2_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[3]$M3 + theme(title = element_blank())
ggsave("../plots/normal_preg/m3_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[4]$M4 + theme(title = element_blank())
ggsave("../plots/normal_preg/m4_interactionplot.png", width = 7, height = 5, dpi = 80)
plots[5]$M5 + theme(title = element_blank())
ggsave("../plots/normal_preg/m5_interactionplot.png", width = 7, height = 5, dpi = 80)


# CEM Heatmap

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

genes <- cem@module[which(cem@module$modules=="M1"),]$genes
m1 <- rep(NA, nrow(m))
names(m1) <- rownames(m)
m1[which(names(m1) %in% genes)] <- mod_cols[["M1"]]

genes <- cem@module[which(cem@module$modules=="M5"),]$genes
m5 <- rep(NA, nrow(m))
names(m5) <- rownames(m)
m5[which(names(m5) %in% genes)] <- mod_cols[["M5"]]

genes <- ltEP$hsymbol
ltEPanno <- rep(NA, nrow(m))
names(ltEPanno) <- rownames(m)
ltEPanno[which(names(ltEPanno) %in% genes)] <- "pink"

genes <- RDP$hsymbol
RDPanno <- rep(NA, nrow(m))
names(RDPanno) <- rownames(m)
RDPanno[which(names(RDPanno) %in% genes)] <- "lightblue"
  
genes <- RAP$hsymbol
RAPanno <- rep(NA, nrow(m))
names(RAPanno) <- rownames(m)
RAPanno[which(names(RAPanno) %in% genes)] <- "lightgreen"
  
rspl <- data.frame(modules = cem@module$modules)
rspl$hsymbol <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M5", "M2", "M4", "M3", "M1", "NC"))

annoRow <- list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5, ltEP = ltEPanno, RDP = RDPanno, RAP = RAPanno)
anno_df <- data.frame(M1 = rownames(m), M2 = rownames(m), M3 = rownames(m), M4 = rownames(m), M5 = rownames(m), ltEP = rownames(m), RDP = rownames(m), RAP = rownames(m))
rha = rowAnnotation(df=anno_df, col=annoRow, show_legend = F)

metadata.left$PhenoNames <- factor(metadata.left$PhenoNames, levels = c("WTnp", "WTpreg", "WTpost"))
metadata.left <- metadata.left[order(metadata.left$PhenoNames),]
m <- m[,match(metadata.left$SampleNumber, colnames(m))]

cspl <- metadata.left$PhenoNames

cha = HeatmapAnnotation(Group = metadata.left[,c("PhenoNames")], col=annoColnp, show_legend = F, show_annotation_name = F)

colors = hcl.colors(30, palette = "Purple-Green")
ht = Heatmap(t(scale(t(m))), show_row_names = F, show_row_dend = F, show_column_names = T, cluster_columns = F,
             top_annotation = cha, right_annotation = rha, name = "expr",
             column_split = cspl, row_split = rspl$modules, cluster_row_slices = F,
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 14),
             column_dend_height=unit(10, "mm"),
             heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14))
)

draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
png("../plots/normal_preg/modules_heatmap.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()
