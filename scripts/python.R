### python

python <- read.table("../metadata/python_genes.csv", sep=",", header = T)
python <- merge(python, t2gh, by.x="symbol", by.y="hsymbol")
rat <- merge(genes_pe_lfc1, t2g, by="ens_gene")
rat <- rat[!duplicated(rat$ens_gene),]
rat <- merge(rat, t2gh, by="ens_gene")
#rat <- rat[which(rat$hsymbol!=""),]
write.table(rat, "../metadata/rat_genes.tsv", sep="\t", row.names = F)

fname = "../plots/final/python_venn.png"
n1 = length(python$ens_gene)
n2 = length(genes_pe_lfc1$ens_gene)

names = c(paste("python, n=", n1, sep=""),
          paste("rat, n=", n2, sep=""))
makeVenn(2, list(python$ens_gene, genes_pe_lfc1$ens_gene),
         names,
         "Venn diagram of python vs rat DEGs",
         fname, c("lightblue", "pink"), dist=0.15, width=950, height=850)

common <- intersect(python$ens_gene, genes_pe_lfc1$ens_gene)
common <- t2gh[which(t2gh$ens_gene %in% common),]
#write.table(common, "../metadata/rat_python_genes.tsv", sep="\t", row.names = F)


pcommon <- python[which(python$symbol %in% common$hsymbol),]
pcommon$log2FC <- -pcommon$log2FC
pcommon <- pcommon[order(pcommon$log2FC),]
pcommon <- pcommon[,c("ens_gene", "symbol", "log2FC", "padj")]
pcommon$comparison <- "fed_starved"


rcm1 <- PEd21_SDd21[which(PEd21_SDd21$ens_gene %in% common$ens_gene),]
rcm1$comparison <- "PEpreg_WTpreg"
rcm2 <- PEpp_SDpp[which(PEpp_SDpp$ens_gene %in% common$ens_gene),]
rcm2$comparison <- "PEpost_WTpost"
rcommon <- rbind(rcm1[,c("ens_gene", "symbol", "log2FoldChange", "padj", "comparison")], rcm2[,c("ens_gene", "symbol", "log2FoldChange", "padj", "comparison")])
colnames(rcommon)[3] <- "log2FC"
rcommon <- merge(rcommon, t2gh, by="ens_gene")
rcommon$symbol <- rcommon$hsymbol
rcommon <- rcommon[,1:5]

allpython <- rbind(pcommon, rcommon)
allpython$sign <- "up"
allpython$sign[which(allpython$log2FC<0)] <- "down"
allpython$sign[which(allpython$padj>0.05)] <- "ns"

p<-ggplot(data=allpython, aes(x=symbol, y=log2FC, fill=sign)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_x_discrete(limits=pcommon$symbol) +
  scale_fill_manual(values=c("blue", "gray", "red")) +
  facet_wrap(~comparison, scales = "free_x")
p

allpython$module <- NA
allpython$module[which(allpython$ens_gene %in% m1$ens_gene)] <- "M1"
allpython$module[which(allpython$ens_gene %in% m2$ens_gene)] <- "M2"
allpython$module[which(allpython$ens_gene %in% m3$ens_gene)] <- "M3"
allpython$module[which(allpython$ens_gene %in% m4$ens_gene)] <- "M4"
allpython$module[which(allpython$ens_gene %in% m5$ens_gene)] <- "M5"

p<-ggplot(data=allpython[which(allpython$comparison=="PEpreg_WTpreg"),], aes(x=symbol, y=log2FC, fill=sign)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_x_discrete(limits=pcommon$symbol) +
  scale_fill_manual(values=c("blue", "red")) +
  facet_wrap(~comparison+module, scales = "free_x")
p

ggVennDiagram(list(python=python$ens_gene, rat=genes_pe_lfc1$ens_gene)) +
  scale_x_continuous(expand = expansion(mult = .2))

which(common$ens_gene %in% m1$ens_gene)
which(common$ens_gene %in% m2$ens_gene)
which(common$ens_gene %in% m3$ens_gene)
which(common$ens_gene %in% m4$ens_gene)
which(common$ens_gene %in% m5$ens_gene)

pythondf <- data.frame(Module=character(), Category=character(), Genes=I(list()), stringsAsFactors=FALSE)
for (i in 1:nrow(ora)) {
  ind <- which(ora$geneID[[i]] %in% common$ens_gene)
  if (length(ind)>0) {
    genes <- ora$geneID[[i]][ind]
    symbols <- rat$hsymbol[which(rat$ens_gene %in% genes)]
    pythondf <- rbind(pythondf, data.frame(Module=ora$Module[i], Category=ora$ID[i], Genes=I(list(c(symbols)))))
  }
}

write.table(pythondf, "../metadata/python_genes_ora.tsv", sep="\t", row.names = F)
