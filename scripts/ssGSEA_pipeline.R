library(ssGSEA2)
library(cmapR)
library(genefilter)

runssGSEA <- function(metadata, abundance, t2gh, filename, path, threshold = 0.5, rerun=FALSE) {
  m <- abundance[, which(colnames(abundance) %in% metadata$SampleNumber)]
  m <- as.data.frame(m)
  m$ens_gene <- rownames(m)
  m <- merge(m, t2gh, by="ens_gene")
  m <- m[!duplicated(m$hsymbol),]
  m <- m[m$hsymbol!="",]
  rownames(m) <- m$hsymbol
  m <- as.matrix(m[,2:(ncol(m)-1)])
  
  cdesc <- data.frame(id=metadata$SampleNumber, type=as.character(metadata$PhenoNames))
  (gct.data <- new("GCT", mat=m, cdesc=cdesc))
  
  write_gct(gct.data, paste(path, filename, ".gct", sep=""), appenddim = FALSE)
  
  if (rerun) {
    res = run_ssGSEA2(paste(path, filename, ".gct", sep=""),
                      output.prefix = filename,
                      gene.set.databases = "../databases/msigdb/h.all.v2023.2.Hs.symbols.gmt",
                      output.directory = paste(path, filename, sep=""),
                      sample.norm.type = "none", 
                      weight = 0.75, 
                      correl.type = "z.score", 
                      statistic = "area.under.RES",
                      output.score.type = "NES", 
                      nperm = 1000, 
                      min.overlap = 25, 
                      extended.output = TRUE, 
                      global.fdr = FALSE,
                      log.file = paste(path, filename, "/run.log", sep=""),
                      spare.cores = 24,
                      par=T)
  }
  
  gct_file <- parse_gctx(paste(path, filename, "/", filename, "-fdr-pvalues.gct", sep=""))
  print(gct_file@mat)
  pthws <- names(which(rowSums(gct_file@mat) < 0.05*ncol(gct_file@mat)))
  
  gct_file <- parse_gctx(paste(path, filename, "/", filename, "-scores.gct", sep=""))
  gct_file@mat
  
  gsea_res <- gct_file@mat[pthws,]
  #print(gsea_res)
  
  gsea_pval <- rowttests(gsea_res, factor(metadata.left$PhenoNames), tstatOnly = FALSE)
  gsea_pval <- gsea_pval[which(gsea_pval$p.value<=0.05),]
  gsea_pval <- gsea_pval[which(abs(gsea_pval$dm)>=threshold),]
  print(gsea_pval)
  
  gct_file <- parse_gctx(paste(path, filename, "/", filename, "-combined.gct", sep=""))
  gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% rownames(gsea_pval)))
  gct_file@mat
  write_gct(gct_file, paste(path, filename, "/", filename, "-combined_filtered.gct", sep=""))
}

metadata.left <- metadata.h[which(grepl("SDd21",metadata.h$Pheno) | grepl("PEd21",metadata.h$Pheno)), ]
metadata.left <- metadata.left[which(metadata.left$Region == "LV"), ]
metadata.left <- metadata.left[-which(metadata.left$SampleNumber %in% c(13,18)),]
filename = "d21"
path = "../degs/ssGSEA/"
runssGSEA(metadata.left, tpm_abd, t2gh, filename, path, threshold = 0.4, rerun=T)

metadata.left <- metadata.h[which(grepl("SDpp",metadata.h$Pheno) | grepl("PEpp",metadata.h$Pheno)), ]
metadata.left <- metadata.left[which(metadata.left$Region == "LV"), ]
filename = "pp"
path = "../degs/ssGSEA/"
runssGSEA(metadata.left, tpm_abd, t2gh, filename, path, threshold = 0.2, rerun=F)


path = "../degs/ssGSEA/pp/signature_gct/"
filename = "HALLMARK_INTERFERON_ALPHA_RESPONSE_n15x94.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% pp$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/pp/signature_filtered/", "HALLMARK_INTERFERON_ALPHA_RESPONSE_n15x6.gct", sep=""))

path = "../degs/ssGSEA/d21/signature_gct/"
filename = "HALLMARK_INTERFERON_ALPHA_RESPONSE_n13x94.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% d21$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/d21/signature_filtered/", "HALLMARK_INTERFERON_ALPHA_RESPONSE_n13x13.gct", sep=""))

path = "../degs/ssGSEA/d21/signature_gct/"
filename = "HALLMARK_ANGIOGENESIS_n13x34.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% d21$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/d21/signature_filtered/", "HALLMARK_ANGIOGENESIS_n13x12_f.gct", sep=""))

path = "../degs/ssGSEA/pp/signature_gct/"
filename = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_n15x192.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% pp$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/pp/signature_filtered/", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_n15x3.gct", sep=""))

path = "../degs/ssGSEA/d21/signature_gct/"
filename = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_n13x192.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% d21$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/d21/signature_filtered/", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_n13x71.gct", sep=""))

path = "../degs/ssGSEA/pp/signature_gct/"
filename = "HALLMARK_MYC_TARGETS_V2_n15x54.gct"
gct_file <- parse_gctx(paste(path, filename, sep=""))
gct_file <- subset_gct(gct_file, rid = which(gct_file@rid %in% pp$hsymbol))
gct_file
write_gct(gct_file, paste("../degs/ssGSEA/pp/signature_filtered/", "HALLMARK_MYC_TARGETS_V2_n15x8.gct", sep=""))
