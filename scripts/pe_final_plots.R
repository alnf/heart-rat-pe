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


### Plot enrichment
get_ora <- function(gmt_path, t2gh, cem) {
  gmt <- read_gmt(gmt_path)
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
    g$grobs[[i]]$children[[2]]$grobs[[2]]$children[[1]]$gp <- gpar(col = colors, lineheight = 0.9)    
    labels <- sub(".*? ", "", labels)
    g$grobs[[i]]$children[[2]]$grobs[[2]]$children[[1]]$label <- labels
  }
  return(as_ggplot(g))
}

#png("../plots/final/enrich_genes_pe_lfc1.png",  width = 10, height = 10, units="in", res=80)
png("../plots/final/enrich_genes_pe_lfc1.png",  width = 10, height = 10, units="in", res=80)
axis_text_color(p, plot_data)
pcols <- c("HALLMARK"="black", "GO Biological Process"="#217364")
lgd = Legend(labels = names(pcols), title = "Pathway", labels_gp = gpar(fontsize = 8), nrow = 1, legend_gp = gpar(fill = pcols))
draw(lgd, x = unit(0.3, "in"), y = unit(0.3, "in"), just = c("left", "bottom"))
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

#draw(us, padding = unit(c(2, 70, 2, 2), "mm"))
#M1 40, 55
#M2 65 190
#M3 25 115
#M4 25 55
#M5 25 15
png(paste("../plots/final/final_upset_genes_pe_lfc1_",modname,".png", sep=""), width = 15, height = 5, units="in", res=100)
draw(us, padding = unit(c(2, 25, 2, 2), "mm"))
dev.off()


oras <- ora[which(ora$Module==modname),]
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
set.seed(42)
cl <- hdbscan(t(ltm), minPts = 2, gen_simplified_tree = F)
oras$Cluster <- as.factor(cl$cluster)
oras$ID <- stringr::str_to_title(oras$ID)
oras$ID <- sub(".*? ", "", oras$ID)
oras$ID <- stringr::str_wrap(oras$ID, 25)


autoplot(pca, data = oras, colour = "Cluster",
         label.show.legend = F) +
  geom_text_repel(aes(label = ID, colour = Cluster), size = 3, max.overlaps = 25, force = 2)
