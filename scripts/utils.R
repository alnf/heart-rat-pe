getDEGs <- function(dds, contrast, t2g, lFC = 1, sval.filter = TRUE, filename, filename_short) {

  res <- results(dds, contrast=contrast)
  res <- lfcShrink(dds, contrast=contrast, type="ashr", lfcThreshold = log2(1.2), alpha=0.05)
  resTable <- data.frame(res)
  resTable$baseMean <- round(resTable$baseMean, 4)
  
  baseMeanPerLvl <- sapply(levels(dds$Group),
                           function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$Group == lvl, drop=F]))   
  baseMeanPerLvl <- baseMeanPerLvl[, contrast[2:3]]
  baseMeanPerLvl <- round(baseMeanPerLvl, 4)

  colnames(baseMeanPerLvl) <- paste("mean", colnames(baseMeanPerLvl), sep="_")
  resTable <- cbind(baseMeanPerLvl, resTable)
  resTable$ens_gene <- rownames(resTable)
  
  resTable <- merge(t2g, resTable, by="ens_gene")
  resTable <- resTable[!duplicated(resTable$ens_gene),]
  resTable$svalue[which(resTable$svalue<0)] <- 0
  resTable <- resTable %>% arrange(desc(abs(log2FoldChange)), padj)

  write.table(resTable, filename, sep="\t", row.names = F)
  
  resTable <- resTable %>% filter(padj <= 0.05)
  if (sval.filter) {
    resTable <- resTable %>% filter(svalue <= 0.005)
  }
  resTable <- resTable %>% filter(!str_detect(description, 'ribosomal RNA'))
  resTable <- resTable %>% filter(symbol!="" | msymbol!="") 
  resTable$symbol[which(resTable$symbol=="")] <- resTable$msymbol[which(resTable$symbol=="")]
  
  mean <- rowMeans(resTable[,6:7])
  if (length(which(mean<20))>0) {
    resTable <- resTable[-which(mean<20),]
  }
  resTable <- resTable %>% filter(abs(log2FoldChange) > lFC)
  print(nrow(resTable))
  write.table(resTable, filename_short, sep="\t", row.names = F)
}

makeHeatmap <- function(dds, contrast, degs, lFC, filename, region) {
  groups <- paste(region, contrast, sep="_")
  ind = which(dds$Group %in% groups)
  
  exprs <- assay(vst(dds[,ind], blind=TRUE))
  colnames(exprs) <- dds$SampleNumber[ind]
  degs.filtered <- degs %>% filter(abs(log2FoldChange) > lFC)
  
  exprs.short <- exprs[degs.filtered$ens_gene,]
  rownames(exprs.short) <- degs.filtered$symbol
  
  color = hue_pal()(2)
  anno.df <- as.data.frame(colData(dds[,ind]))
  rownames(anno.df) <- anno.df$SampleNumber
  anno.df$Group <- factor(anno.df$Group, levels=groups)
  print(anno.df$Group)
  
  names(color) <- levels(factor(anno.df$Group))
  annoCol <- list(Group = color)

  pheatmap(exprs.short, annotation = anno.df[,c("Group"), drop=F], annotation_colors = annoCol,
           filename = filename, scale = "row", width = 8, height = 11,
           main = paste(contrast[1], "vs", contrast[2], ": n =", nrow(degs.filtered)))
    
}

makeHeatmapFC <- function(fc, genes, contrast, filename, scale = "row") {
  rownames(fc) <- genes$symbol

  color = hue_pal()(4)
  anno.df <- data.frame(region = colnames(fc))
  rownames(anno.df) <- anno.df$region
  anno.df$region <- factor(anno.df$region, levels=colnames(fc))
  names(color) <- levels(factor(anno.df$region))
  annoCol <- list(region = color)
  
  
  pheatmap(fc, annotation = anno.df, annotation_colors = annoCol,
           filename = filename, scale = scale, width = 8, height = 11,
           main = paste(contrast[1], "vs", contrast[2], " all regions: n =", nrow(fc)))
  
}



makeVenn <- function(n, genesList, names, title, filename, colors, dist = 0.02, width=520, height=480){
  venn.diagram(
    x = genesList,
    category.names = names,
    filename = filename,
    output=TRUE,
    disable.logging = TRUE,
    main = title,
    
    # Output features
    imagetype="png" ,
    height = height , 
    width = width , 
    resolution = 300,
    compression = "lzw",
    margin = 0.05, 

    # Circles
    lwd = 2,
    lty = 'blank',
    fill = colors,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    ext.text = FALSE,
    #ext.line.lwd = 2,
    #ext.dist = -0.15,
    #ext.length = 0.9,
    #ext.pos = -4,
    
    # Set names
    main.cex = 0.4,
    main.fontface = "bold",
    cat.cex = 0.5,
    #cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = rep(0, n),
    cat.dist = rep(dist, n)
  )
}

genePlot <- function(dds, gene, intgroup, groups, comparisons, filename, title){
  d <- plotCounts(dds, gene=which(rownames(dds)==gene), intgroup=intgroup, 
                  returnData=TRUE, normalized = TRUE, transform = TRUE)
  
  d <- d[which(d$Group %in% groups),]
  d$count <- log2(d$count)
  ggboxplot(data=d, x="Group", y="count", color="Group", palette = "jco", add = "jitter", title = title) +
    stat_compare_means(method = "t.test", comparisons = comparisons)
  ggsave(filename, width = 7, height = 5, dpi=100)
}

genePlotRegionsBar <- function(dds, gene, intgroup, comparisons, filename, title){
  d <- plotCounts(dds, gene=which(rownames(dds)==gene), intgroup=intgroup, 
                  returnData=TRUE, normalized = TRUE, transform = TRUE)
  
  d$count <- log2(d$count)
  d$Region <- factor(sapply(strsplit(as.character(d$Group),"_"), `[`, 1))
  d$Pheno <- factor(sapply(strsplit(as.character(d$Group),"_"), `[`, 2))
  
  ggboxplot(data=d, x="Pheno", y="count", color="Pheno", palette = "jco", add = "jitter", title = title) +
    facet_wrap(~Region, nrow = 2, ncol = 2) +
    stat_compare_means(method = "t.test", comparisons = comparisons) +
    theme_bw()

  ggsave(filename, width = 10, height = 7, dpi=110)
}

genePlotRegionsLine <- function(dds, gene, intgroup, comparisons, filename, title){
  d <- plotCounts(dds, gene=which(rownames(dds)==gene), intgroup=intgroup, 
                  returnData=TRUE, normalized = TRUE, transform = TRUE)
  
  d$count <- log2(d$count)
  d$Region <- factor(sapply(strsplit(as.character(d$Group),"_"), `[`, 1))
  d$Pheno <- factor(sapply(strsplit(as.character(d$Group),"_"), `[`, 2))
  
  ggboxplot(data=d, x="Pheno", y="count", color="Pheno", palette = "jco", add = "jitter", title = title) +
    facet_wrap(~Region, nrow = 4) +
    stat_compare_means(method = "t.test", comparisons = comparisons) +
    theme_bw()
  
  d.sum <- d %>%
    group_by(Region, Pheno) %>%
    summarise(
      sd = sd(count, na.rm = TRUE),
      mcount = mean(count)
    )
  
  ggplot(data=d.sum, aes(x=Pheno, y=mcount, group=Region, color=Region)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mcount-sd, ymax=mcount+sd), width=.2,
                  position=position_dodge(0.2))  
  
  ggsave(filename, width = 8, height = 6, dpi=100)
}


regionPlot <- function(dds, region, ens_genes, filename, width=10, height=7){
  ind <- which(grepl(region, dds$Group))
  anno.df <- as.data.frame(colData(dds[,ind]))
  rownames(anno.df) <- anno.df$SampleNumber
  
  exprs <- assay(vst(dds[,ind], blind=TRUE))
  exprs <- exprs[ens_genes,]
  colnames(exprs) <- anno.df$SampleNumber
  cor <- cor(exprs)
  color = hue_pal()(length(levels(factor(dds[,ind]$Group))))
  names(color) <- levels(factor(anno.df$Group))
  annoCol <- list(Group = color)
  pheatmap(cor, annotation = anno.df[,c("Group"), drop=F], annotation_colors = annoCol,
           filename = filename, width = width, height = height,
           main = paste(region, ", ngenes=", length(ens_genes), sep=""))    
}

MAplot <- function(resDegs, mean_limit, fc_limit, nudge, title){
  ggplot(degs, aes(x=log2(baseMean+1), y=log2FoldChange, color=degs, label=hsymbol)) +
    geom_point(alpha=0.7) +
    geom_hline(yintercept=0, color = "gray", size=1) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    geom_text_repel(data          = subset(degs, (log2(baseMean+1) > mean_limit) & !is.na(degs) & abs(log2FoldChange) > fc_limit),
                    nudge_y       = nudge,
                    nudge_x       = nudge,
                    size          = 6,
                    box.padding   = 0.5,
                    point.padding = 0.5,
                    force         = 1,
                    segment.color = "grey50",
                    direction     = "both",
                    max.overlaps = 20,
                    show.legend = F) +
    labs(title=title) +
    theme(plot.title = element_text(vjust=-7, hjust=1, size=16),
          axis.title = element_text(size=14),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=10))
}


compHeatmap <- function(pheno, region, order=NA, metadata, lcounts, glist, csplit=T, dlists, cols, annoCol, clegend=F, clcol=T){
  meta <- metadata[which(metadata$Region == region), ]
  meta <- meta[which(meta$Pheno %in% pheno), ]
  m <- lcounts[glist$ens_gene, which(colnames(lcounts) %in% meta$SampleNumber)]

  annoRow <- list()
  for (i in 1:length(dlists)) {
    genes <- dlists[[i]]$ens_gene
    anno <- rep(NA, nrow(m))
    names(anno) <- rownames(m)
    anno[which(names(anno) %in% genes)] <- cols[[names(dlists)[i]]]
    anno <- list(anno) 
    names(anno) <- names(dlists)[i]
    annoRow <- append(annoRow, anno)
  }
  
  if (!is.null(order)) {
    meta$PhenoNames <- factor(meta$PhenoNames, levels = order)
    meta <- meta[order(meta$PhenoNames),]
    m <- m[,match(meta$SampleNumber, colnames(m))]
  }
  
  if (csplit) {
    cspl <- meta$PhenoNames
  } else {
    cspl = NULL
  }
  
  anno_df = data.frame(matrix(NA, nrow = nrow(m), ncol = length(names(dlists))))
  for (i in 1:ncol(anno_df)) {
    anno_df[,i] <- rownames(m)
  }
  colnames(anno_df) <- names(dlists)
  
  rha = rowAnnotation(df=anno_df, col=annoRow, show_legend = F)
  cha = HeatmapAnnotation(Group = meta[,c("PhenoNames")], col=annoCol, show_legend = clegend, show_annotation_name = F)
  
  
  ht = Heatmap(t(scale(t(m))), show_row_names = F, show_row_dend = F, show_column_names = T, cluster_columns = clcol,
               top_annotation = cha, right_annotation = rha, name = "expr",
               column_split = cspl, row_split = NULL, cluster_row_slices = T, cluster_column_slices = T,
               row_names_gp = gpar(fontsize = 14),
               column_names_gp = gpar(fontsize = 14),
               column_dend_height=unit(10, "mm"),
               heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14))
  )
  return(ht)
}
