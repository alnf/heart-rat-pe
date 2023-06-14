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

makeHeatmap <- function() {
  ind = which(dds$Group %in% contrast[2:3])
  exprs <- assay(vst(dds[,ind], blind=TRUE))
  colnames(exprs) <- dds$SampleNumber[ind]
  
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
  write.table(resTable, paste("../degs/WT/", region, "/degs_", fname, "_short.tsv", sep=""), sep="\t", row.names = F)
  
  resTable <- resTable %>% filter(abs(log2FoldChange) > lFCvis)
  exprs.short <- exprs[resTable$ens_gene,]
  cor <- cor(exprs.short)
  color = hue_pal()(2)
  anno.df <- as.data.frame(colData(dds[,ind]))
  rownames(anno.df) <- anno.df$SampleNumber
  names(color) <- levels(factor(anno.df$Group))
  annoCol <- list(Group = color)
  pheatmap(cor, annotation = anno.df[,c("Group"), drop=F], annotation_colors = annoCol,
           filename = paste("../plots/WT/", region, "/heatmap_vst_cor_", fname, ".png", sep=""),
           width = width, height = height)    
  rownames(exprs.short) <- resTable$symbol
  pheatmap(exprs.short, annotation = anno.df[,c("Group"), drop=F], annotation_colors = annoCol,
           filename = paste("../plots/WT/", region, "/heatmap_vst_", fname, ".png", sep=""),
           scale = "row", width = 8, height = 11,
           main = paste(contrast[2], "vs", contrast[3], ": n =", nrow(resTable)))
    
}


makeVenn <- function(n, genesList, names, title, filename, colors, dist = 0.02){
  venn.diagram(
    x = genesList,
    category.names = names,
    filename = filename,
    output=TRUE,
    disable.logging = TRUE,
    main = title,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
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
    main.cex = 0.6,
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
  ggsave(filename, width = 7, height = 5)
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

