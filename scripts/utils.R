getDEGs <- function(dds, contrast, t2g, width, height, lFC = 1, lFCvis = 2) {

  #res <- results(dds, contrast=contrast, lfcThreshold = 0.5, altHypothesis = "greaterAbs", alpha=0.05)
  res <- results(dds, contrast=contrast)
  res <- lfcShrink(dds, contrast=contrast, type="ashr", lfcThreshold = log2(1.2),
                   alpha=0.05)
  summary(res)

  resTable <- data.frame(res)
  resTable <- resTable %>% filter(padj <= 0.05)

  baseMeanPerLvl <- sapply(levels(dds$Group),
                           function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$Group == lvl, drop=F]))   
  baseMeanPerLvl <- baseMeanPerLvl[rownames(resTable), contrast[2:3]]
  baseMeanPerLvl <- round(baseMeanPerLvl, 4)

  colnames(baseMeanPerLvl) <- paste("mean", colnames(baseMeanPerLvl), sep="_")
  resTable <- cbind(baseMeanPerLvl, resTable)
  resTable$ens_gene <- rownames(resTable)
  
  resTable <- merge(t2g, resTable, by="ens_gene")
  resTable <- resTable[!duplicated(resTable$ens_gene),]
  resTable <- resTable %>% arrange(desc(abs(log2FoldChange)))
  resTable$svalue[which(resTable$svalue<0)] <- 0
  
  fname = paste(contrast[2:3], collapse="_")
  region = strsplit(contrast[2], "_")[[1]][1]
  write.table(resTable, paste("../degs/WT/", region, "/degs_", fname, ".tsv", sep=""), sep="\t", row.names = F)

  ind = which(dds$Group %in% contrast[2:3])
  exprs <- assay(vst(dds[,ind], blind=TRUE))
  colnames(exprs) <- dds$SampleNumber[ind]

  resTable <- resTable %>% filter(svalue <= 0.005)  
  ttt <- resTable
  resTable <- ttt
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
