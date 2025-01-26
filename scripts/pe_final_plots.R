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

metadata.left <- metadata.h[which(metadata.h$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$Pheno == "np"), ]

t2gh <- biomaRt::getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), mart = ensembl)
t2gh <- dplyr::rename(t2gh, ens_gene = ensembl_gene_id, hsymbol = hsapiens_homolog_associated_gene_name)

### All genes
types = c("_short.tsv", ".tsv")
type = types[2]

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

#fviz_eig(pca, addlabels = TRUE, ncp=10)
#lbls <- PEd21_SDd21$symbol
#names(lbls) <- PEd21_SDd21$ens_gene
#fviz_contrib(pca, choice = "var", axes = 1, top = 20) + scale_x_discrete(labels = lbls)
#fviz_contrib(pca, choice = "var", axes = 2, top = 20) + scale_x_discrete(labels = lbls)

## Plotting Venns

names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/final/final_venn_", fname, ".png", sep="")
n1 = length(PEd21_SDd21$ens_gene[which(abs(PEd21_SDd21$log2FoldChange)>1)])
n2 = length(PEpp_SDpp$ens_gene[which(abs(PEpp_SDpp$log2FoldChange)>1)])
n3 = length(PEd21_PEpp$ens_gene[which(abs(PEd21_PEpp$log2FoldChange)>1)])
n4 = length(SDd21_SDpp$ens_gene[which(abs(SDd21_SDpp$log2FoldChange)>1)])

names = c(paste("preg, n=", n1, sep=""),
          paste("post, n=", n2, sep=""),
          paste("deltaPE, n=", n3, sep=""),
          paste("deltaWT, n=", n4, sep=""))
makeVenn(4, list(PEd21_SDd21$ens_gene[which(abs(PEd21_SDd21$log2FoldChange)>1)], PEpp_SDpp$ens_gene[which(abs(PEpp_SDpp$log2FoldChange)>1)],
                 PEd21_PEpp$ens_gene[which(abs(PEd21_PEpp$log2FoldChange)>1)], SDd21_SDpp$ens_gene[which(abs(SDd21_SDpp$log2FoldChange)>1)]),
         names,
         paste("Venn diagram of DEGs in 4 groups of comparison,\n lFC=1, ngenes=", nrow(genes_pe_lfc1), sep=""),
         fname, c(mcols["preg"], mcols["post"], mcols["deltaPE"], mcols["deltaWT"]), dist=0.15, width=650, height=550)

## Plotting heatmaps, ora

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
m1n <- nrow(m1)
m2n <- nrow(m2)
m3n <- nrow(m3)
m4n <- nrow(m4)
m5n <- nrow(m5)


#dlists <- list(preg=PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),], post=PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),],
#               deltaPE=PEd21_PEpp[which(abs(PEd21_PEpp$log2FoldChange)>1),], deltaWT=SDd21_SDpp[which(abs(SDd21_SDpp$log2FoldChange)>1),])
dlists <- list(preg=PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),], post=PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),])
cols <- list(preg=mcols["preg"], post=mcols["post"], deltaPE=mcols["deltaPE"], deltaWT=mcols["deltaWT"])
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
pheno <- order
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=F, csplit=T, dlists, cols, acols, clegend=F, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_lfc1_2vars.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()

# With clusterting
order <- c("PEpreg", "PEpost", "WTpreg", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, cc=T, csplit=T, dlists, cols, acols, clegend=T, rsplit=rspl$modules, modcols)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
png("../plots/final/final_heatmap_genes_pe_lfc1_dd.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(5, 2, 2, 3), "mm"), heatmap_legend_side = "bottom")
dev.off()


### Plot enrichment
get_ora <- function(gmt_path, t2gh, cem) {
  gmt <- CEMiTool::read_gmt(gmt_path)
  gmt <- merge(gmt, t2gh, by.x="gene", by.y="hsymbol")
  gmt <- gmt[,c(2:3)]
  colnames(gmt)[2] <- "gene"
  cem <- mod_ora(cem, gmt)
  ora <- cem@ora
  return(ora)
}

ora_hm <- get_ora("../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt", t2gh, cem)
ora_bp <- get_ora("../databases/msigdb/c5.go.bp.v2023.2.Hs.symbols.gmt", t2gh, cem)
ora <- rbind(ora_hm, ora_bp)

ora$pathway <- sapply(strsplit(ora$ID,"_"), `[`, 1)
#ora$ID <- sub(".*?_", "", ora$ID)
ora$ID <- gsub("_", " ", ora$ID)
ora$geneRatio <- sapply(ora$GeneRatio, function(x) eval(parse(text=x)))

ora <- ora[which(ora$Count>5),]
ora %>% group_by(Module) %>% slice_max(order_by = geneRatio, n = 10) -> ora
ora$col <- ifelse(ora$pathway == "HALLMARK", "black", "#217364")

strip <- strip_themed(text_y = elem_list_text(color = rep("white", 5), size=rep(14, 5), face=rep("bold", 5)),
                      background_y = elem_list_rect(fill = mcols[9:13]))

plot_data <- ora %>%
  arrange(Module, geneRatio) %>%
  mutate(rank = row_number())

plot_data$ID <- stringr::str_to_title(plot_data$ID)
plot_data$y <- 1:nrow(plot_data)
plot_data$y <- factor(plot_data$y)

idnames <- as.character(plot_data$ID)
names(idnames) <- plot_data$y

plot_data$pathway <- factor(plot_data$pathway)

p <- ggplot(plot_data, aes(geneRatio, y)) + 
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
      scale_size_continuous(range=c(1, 7)) +
      geom_segment(aes(xend=0, yend = y)) +
      theme_minimal() +
      xlab("Gene Ratio") +
      ylab(NULL) + 
      facet_grid2(Module~., scales="free_y", strip = strip, drop=T, axes = "margins", space = "free_y") +
      scale_y_discrete(labels = idnames) +
      theme(axis.text.y = element_text(size = 11))
p
axis_text_color <- function(plot, plot_data) {
  g <- ggplotGrob(plot)
  yaxis_grobs <- which(grepl("axis-l", g$layout$name))
  for (i in yaxis_grobs) {
    labels <- g$grobs[[i]]$children[[2]]$grobs[[2]]$children[[1]]$label
    colors <- plot_data[which(plot_data$ID %in% labels), c("ID", "col")]
    colors <- colors[!duplicated(colors$ID),]
    colors <- colors[match(labels, colors$ID),]
    colors <- colors$col
    names(colors) <- labels
    #g$grobs[[i]]$children[[2]]$grobs[[2]]$children[[1]]$gp <- gpar(col = colors, lineheight = 0.9)    
    labels <- sub(".*? ", "", labels)
    g$grobs[[i]]$children[[2]]$grobs[[2]]$children[[1]]$label <- labels
  }
  return(as_ggplot(g))
}

png("../plots/final/enrich_genes_pe_lfc1.png",  width = 10, height = 10, units="in", res=80)
axis_text_color(p, plot_data)
#pcols <- c("HALLMARK"="black", "GO Biological Process"="#217364")
#lgd = Legend(labels = names(pcols), title = "Pathway", labels_gp = gpar(fontsize = 8), nrow = 1, legend_gp = gpar(fill = pcols))
#draw(lgd, x = unit(0.3, "in"), y = unit(0.3, "in"), just = c("left", "bottom"))
dev.off()


## Plot upset
ora$geneID <- sapply(strsplit(ora$geneID,"/"), `[`, )
names(ora$geneID) <- ora$ID

modname = "M1"

degree = 1
ups <- make_comb_mat(ora[which(ora$Module==modname),]$geneID, mode = "intersect")
#ups <- ups[comb_size(ups) >= 1 & comb_degree(ups) == degree]

us <- UpSet(ups, pt_size = unit(4, "mm"), lwd = 2,
            comb_col = brewer.pal(comb_degree(ups)[1], "Paired")[comb_degree(ups)],
            top_annotation = upset_top_annotation(ups, add_numbers = TRUE),
            right_annotation = upset_right_annotation(ups, add_numbers = TRUE),
            row_names_gp = grid::gpar(fontsize = 9),
            comb_order = order(comb_size(ups)))

#M1 40, 55
#M2 65 190
#M3 25 115
#M4 25 55
#M5 25 15
png(paste("../plots/final/final_upset_genes_pe_lfc1_",modname,".png", sep=""), width = 15, height = 5, units="in", res=100)
draw(us, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()

## Plot PCA of enrichment
modname = "M5"

oras <- ora[which(ora$Module==modname),]
ltm <- list_to_matrix(oras$geneID)
pca <- prcomp(t(ltm))
rownames(oras) <- oras$ID
set.seed(42)
cl <- hdbscan(t(ltm), minPts = 2, gen_simplified_tree = F)
oras$Cluster <- as.factor(cl$cluster)
oras$ID <- stringr::str_to_title(oras$ID)
oras$ID <- sub(".*? ", "", oras$ID)
oras$ID <- stringr::str_wrap(oras$ID, 25)

rcolors <- c("HALLMARK"="black", "GOBP"="#217364", "0"="#FF0A54", "1"="#4E3D42", "2"="#0A78C2", "3"="#17BEBB")
oras$cluscols <- rcolors[as.character(oras$Cluster)]
names(oras$col) <- oras$pathway
rcols <- list(cluster = oras$cluscols, pathway = oras$col)

autoplot(pca, data = oras, colour = "Cluster",
         label.show.legend = F) +
  geom_text_repel(aes(label = ID, colour = Cluster), size = 3.5, max.overlaps = 25, force = 2.5,
                  direction = "both", max.time = 5, max.iter = 20000, nudge_y = -0.2, segment.color = "gray", segment.alpha=0.5,
                  show.legend = FALSE) +
  scale_colour_manual(name = "Cluster", values = oras$cluscols)

ggsave(paste("../plots/final/pca_enrich_pe_lfc1_",modname,".png", sep=""), width = 8, height = 6, scale = 0.8, dpi = 110)


## Plot heatmap of enrichment

colors = structure(c("#494C6F", "#D80032"), names = c("0", "1"))
rha = rowAnnotation(cluster = oras$Cluster,
                    #pathway = oras$pathway,
                    show_legend = F, col=rcols)


get_fc_heatmap <- function(gnames, degs_list) {
  ltm_df <- data.frame(row.names = rownames(ltm), ens_gene = rownames(ltm))
  for (i in 1:length(degs_list)) {
    ltm_df <- merge(ltm_df, degs_list[[i]][, c("ens_gene", "log2FoldChange", "padj")], by="ens_gene", all.x=T)
    ltm_df[which(ltm_df[,ncol(ltm_df)]>0.05),ncol(ltm_df)-1] <- NA
    ltm_df <- ltm_df[,-c(ncol(ltm_df))]
  }
  rownames(ltm_df) <- ltm_df$ens_gene
  ltm_df <- ltm_df[,-c(1)]
  colnames(ltm_df) <- names(degs_list)
  pch = rep("x", nrow(ltm))
  return(ltm_df)
}


degs_list <- list(PEd21_SDd21, PEpp_SDpp)
names(degs_list) <- c("preg", "post")
ltm_df1 <- get_fc_heatmap(rownames(lts), degs_list)
pch1 <- as.matrix(ifelse(abs(ltm_df1)>1, "‧", NA))
pch1col <- as.matrix(ifelse(ltm_df1>1, "black", NA))
pch1col[which(ltm_df1<(-1))] <- "white"

  
degs_list <- list(PEd21_PEpp, SDd21_SDpp)
names(degs_list) <- c("deltaPE", "deltaWT")
ltm_df2 <- get_fc_heatmap(rownames(lts), degs_list)
pch2 <- as.matrix(ifelse(abs(ltm_df2)>1, "‧", NA))
pch2col <- as.matrix(ifelse(ltm_df2>1, "black", NA))
pch2col[which(ltm_df2<(-1))] <- "white"
  
min_ltm <- min(min(ltm_df1, na.rm=T), min(ltm_df2, na.rm=T))
max_ltm <- max(max(ltm_df1, na.rm=T), max(ltm_df2, na.rm=T))

if (min_ltm>0) {
  col_fun = colorRamp2(c(0, max_ltm), c("white", "red"))
} else if (max_ltm<0) {
  col_fun = colorRamp2(c(min_ltm, 0), c("blue","white"))
} else {
  col_fun = colorRamp2(c(min_ltm, 0, max_ltm), c("blue", "white", "red"))
}


lgd1 = Legend(col_fun = col_fun, title = "log2FC", direction="horizontal")
lgd2 = Legend(labels=c("logFC>1", "logFC<-1", "not significant"), direction="horizontal",
              pch="‧", type = "points",
              legend_gp = gpar(col = c("black", "white", "black")),
              background = c("red","blue", "black"))

gnames <- data.frame(ens_gene=rownames(ltm))
clabs <- merge(gnames, t2g, by="ens_gene")
clabs <- clabs[!duplicated(clabs$ens_gene),]
clabs[which(clabs$symbol==""),]$symbol <- clabs[which(clabs$symbol==""),]$msymbol

ind <- c()
ind <- c(ind, which(clabs$symbol=="Dpp4"))
ind <- c(ind, which(clabs$symbol=="Mmp12"))
ind <- c(ind, which(clabs$symbol=="Trem2"))
ind <- c(ind, which(clabs$symbol=="Postn"))
ind <- c(ind, which(clabs$symbol=="Cd1d1"))
ind <- c(ind, which(clabs$symbol=="Gbp1"))
ind <- c(ind, which(clabs$symbol=="Vegfd"))
ind <- c(ind, which(clabs$symbol=="Rrm2"))
ind <- c(ind, which(clabs$symbol=="Mx1"))


clabs$symbol[ind] <- paste(clabs$symbol[ind], "────────", sep=" ")
fontfaces = rep("plain", nrow(ltm))
fontfaces[ind] <- "bold"


cha1 = HeatmapAnnotation(preg=anno_simple(ltm_df1$preg, pch = pch1[,1], col=col_fun, na_col = "black", pt_gp=gpar(col=pch1col[,1])),
                         post=anno_simple(ltm_df1$post, pch = pch1[,2], col=col_fun, na_col = "black", pt_gp=gpar(col=pch1col[,2])),
                         show_annotation_name = T, annotation_name_side = "left",
                         genes=anno_text(clabs$symbol, gp = gpar(fontface = fontfaces, fontsize=11)),
                         gap = unit(2, "points"))

cha2 = HeatmapAnnotation(deltaWT=anno_simple(ltm_df2$deltaWT, pch = pch2[,2], col=col_fun, na_col = "black", pt_gp=gpar(col=pch2col[,2])),
                         deltaPE=anno_simple(ltm_df2$deltaPE, pch = pch2[,1], col=col_fun, na_col = "black", pt_gp=gpar(col=pch2col[,1])),
                         show_annotation_name = T, annotation_name_side = "left",
                         gap = unit(2, "points"))

ht = Heatmap(t(ltm), clustering_distance_rows = "binary", right_annotation = rha, bottom_annotation = cha1, top_annotation = cha2,
             row_split = oras$Cluster,
             name = "gene", rect_gp = gpar(col = "#725A62", lwd = 1), show_heatmap_legend = F,
             col = colors, row_names_max_width = unit(3000, "mm"),
             row_labels = stringr::str_wrap(oras$ID, 30),
             column_labels = clabs$symbol,
             column_names_rot = 90,
             column_names_gp = gpar(fontsize = 11),
             show_column_names = F)

draw(ht, padding = unit(c(7, 2, 2, 10), "mm"), annotation_legend_list = list(lgd1, lgd2), annotation_legend_side = "bottom")

if (modname == "M1") {
  w = 25
  h = 7.8
}
if (modname == "M2") {
  w = 11
  h = 9
}
if (modname == "M3") {
  w = 13
  h = 7
}
if (modname == "M4") {
  w = 7
  h = 7
}
if (modname == "M5") {
  w = 7
  h = 6
}
png(paste("../plots/final/heatmap_enrich_pe_lfc1_",modname,".png", sep=""), width = w, height = h, units="in", res=100)
draw(ht, padding = unit(c(7, 2, 2, 10), "mm"), annotation_legend_list = list(lgd1, lgd2), annotation_legend_side = "bottom")
dev.off()

#sclist <- read.table("../metadata/sc_list.tsv", sep="\t", header = T)
#sclist <- sclist[!duplicated(sclist$hsymbol),,drop=F]
#sclist <- merge(sclist, t2gh, by="hsymbol")
#sclist <- sclist[!duplicated(sclist$ens_gene),]


# MA plots
imp <- PEd21_SDd21[which(abs(PEd21_SDd21$log2FoldChange)>1),]
degs <- read.table("../degs/LV/degs_LV_PEd21_LV_SDd21.tsv", sep="\t", header = T, stringsAsFactors = F)

imp <- PEpp_SDpp[which(abs(PEpp_SDpp$log2FoldChange)>1),]
degs <- read.table("../degs/LV/degs_LV_PEpp_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)

imp <- PEd21_PEpp[which(abs(PEd21_PEpp$log2FoldChange)>1),]
degs <- read.table("../degs/LV/degs_LV_PEd21_LV_PEpp.tsv", sep="\t", header = T, stringsAsFactors = F)

imp <- SDd21_SDpp[which(abs(SDd21_SDpp$log2FoldChange)>1),]
degs <- read.table("../degs/LV/degs_LV_SDd21_LV_SDpp.tsv", sep="\t", header = T, stringsAsFactors = F)

imp$degs <- NA
imp$degs[which(imp$log2FoldChange>0)] <- "up"
imp$degs[which(imp$log2FoldChange<0)] <- "down"

degs <- merge(degs, imp[c("ens_gene", "degs")], by="ens_gene", all.x=T)
degs <- degs[which(!is.na(degs$log2FoldChange)),]
degs <- degs[-which(degs$ens_gene==agt$ens_gene),]
degs$degs <- factor(degs$degs, levels = c("up", "down"))
degs <- degs %>% filter(!str_detect(description, 'ribosomal RNA'))

MAplot(degs, 7.5, 1, 0.5, paste("PEpreg vs WTpreg, |logFC>1|, ngenes =", nrow(imp)), 15, length(which(imp$degs=="up")), length(which(imp$degs=="down")))
ggsave("../plots/final/ma_preg.png", width = 10, height = 8, dpi = 100)

MAplot(degs, 7.5, 1, 0.5, paste("PEpost vs WTpost, |logFC>1|, ngenes =", nrow(imp)), 15, length(which(imp$degs=="up")), length(which(imp$degs=="down")))
ggsave("../plots/final/ma_post.png", width = 10, height = 8, dpi = 100)

MAplot(degs, 7.5, 1, 0.5, paste("PEpreg vs PEpost, |logFC>1|, ngenes =", nrow(imp)), 15, length(which(imp$degs=="up")), length(which(imp$degs=="down")))
ggsave("../plots/final/ma_deltaPE.png", width = 10, height = 8, dpi = 100)

MAplot(degs, 7.5, 1, 0.5, paste("WTpreg vs WTpost, |logFC>1|, ngenes =", nrow(imp)), 15, length(which(imp$degs=="up")), length(which(imp$degs=="down")))
ggsave("../plots/final/ma_deltaSD.png", width = 10, height = 8, dpi = 100)


## Average expression per module
m <- lcounts[genes_pe_lfc1$ens_gene, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
mods <- cem@module

grouped_columns <- split(seq_along(colnames(m)), factor(metadata.left$PhenoNames))
row_means_by_group <- lapply(grouped_columns, function(cols) rowMeans(m[, cols]))
result <- do.call(cbind, row_means_by_group)
colnames(result) <- names(grouped_columns)
result <- as.data.frame(result)
result$genes <- rownames(result)
result <- merge(result, mods, by="genes")
result <- melt(result, id=c("genes", "modules"), variable.name="group", value.name="expression")

ggplot(result, aes(x=group, y=expression, fill=group, group=genes)) +
  geom_boxplot(aes(group=group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.2) +
  geom_line(alpha=0.5, color="red")+
  facet_wrap(~modules) +
  scale_fill_manual(values = mcols)


## FGSEA




## More ontologies
library(mulea)
library(ExperimentHub)

m1h <- merge(m1, t2gh, by="ens_gene")
m1h <- m1h[!duplicated(m1h$ens_gene),]
colnames(m1h) <- c("rat_ens_gene", "human_symbol")
write.table(m1h, "../metadata/m1_genes.tsv", sep="\t", row.names = F)

which(m1h$hsymbol %in% sclist$hsymbol)

m2h <- merge(m2, t2gh, by="ens_gene")
m2h <- m2h[!duplicated(m2h$ens_gene),]
colnames(m2h) <- c("rat_ens_gene", "human_symbol")
write.table(m2h, "../metadata/m2_genes.tsv", sep="\t", row.names = F)


which(m2h$hsymbol %in% sclist$hsymbol)


tf_ontology <- mulea::read_gmt("../databases/mulea/GO_BP_Rattus_norvegicus_EnsemblID.gmt")
target_set <- m1$ens_gene
background_set  <- rownames(lcounts)

ora_model <- ora(gmt = tf_ontology, 
                 element_names = target_set, 
                 background_element_names = background_set, 
                 p_value_adjustment_method = "eFDR", 
                 number_of_permutations = 10000,
                 nthreads = 24, 
                 random_seed = 42) 

ora_results <- run_test(ora_model)

ora_results %>%
  arrange(eFDR) %>% 
  filter(eFDR < 0.05) -> ora_results

ora_reshaped_results <- reshape_results(model = ora_model, 
                                        model_results = ora_results, 
                                        p_value_type_colname = "eFDR")


plot_lollipop(reshaped_results = ora_reshaped_results,
              # Column containing the names we wish to plot
              ontology_id_colname = "ontology_id",
              # Upper threshold for the value indicating the significance
              p_value_max_threshold = 0.05,
              # Column that indicates the significance values
              p_value_type_colname = "eFDR")



gmt <- read_gmt("../databases/mulea/GO_BP_Rattus_norvegicus_EnsemblID.gmt")


library(msigdbr)

all_gene_sets = msigdbr(species = "Rattus norvegicus")
head(all_gene_sets)
factor(all_gene_sets$gs_cat)


print(msigdbr_collections(), n=23)
msigdbr_collections()


## GMT from Mulea

gmt <- read_gmt("../databases/mulea/GO_BP_Rattus_norvegicus_EnsemblID.gmt")
gmt$list_of_values[1]
gmt <- gmt[,2:3]
colnames(gmt) <- c("term", "gene")
gmt %>% unnest(gene) -> gmt

gmt$gene <- sapply(strsplit(gmt$gene,"\\."), function(x) x[1])
cem <- mod_ora(cem, gmt)

ora <- cem@ora
ora$geneRatio <- sapply(ora$GeneRatio, function(x) eval(parse(text=x)))
ora <- ora[which(ora$Count>5),]
ora %>% group_by(Module) %>% slice_max(order_by = geneRatio, n = 10) -> ora
