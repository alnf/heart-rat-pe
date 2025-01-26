library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

impanel <- c("ENSRNOG00000007468", #IFNG
             "ENSRNOG00000004649", #IL-1B
             "ENSRNOG00000007624", #IL-4
             "ENSRNOG00000008111", #IL-5
             "ENSRNOG00000010278", #IL-6
             "ENSRNOG00000004647", #IL-10
             "ENSRNOG00000007652", #IL-13
             "ENSRNOG00000002802", #CXCL-1
             "ENSRNOG00000070745"  #TNF
             )

impanel <- c("ENSRNOG00000010208", #TIMP1
             "ENSRNOG00000007159", #CCL2
             "ENSRNOG00000056219", #OLR1
             "ENSRNOG00000010665", #CCR7
             "ENSRNOG00000065741", #IL7R
             "ENSRNOG00000012779", #MSR1
             "ENSRNOG00000003622", #CYBB
             "ENSRNOG00000000239", #CCL7
             "ENSRNOG00000014361", #EDN1
             "ENSRNOG00000004273", #IFITM1
             "ENSRNOG00000046254", #ADGRE1
             "ENSRNOG00000004500", #MYC
             "ENSRNOG00000017819", #CD14
             "ENSRNOG00000008759", #Csf3r
             "ENSRNOG00000009211", #C3ar1
             "ENSRNOG00000014320"  #Inhba
)


mleft <- metadata.left[,c("SampleNumber", "PhenoNames")]


imexp <- as.data.frame(lcounts[rownames(lcounts) %in% impanel,])
imexp$ens_gene <- rownames(imexp)
imexp <- melt(imexp, id="ens_gene")
colnames(imexp)[2:3] <- c("SampleNumber", "expr")
imexp <- merge(imexp, mleft, by="SampleNumber")

t2gsimple <- t2g[,c("ens_gene", "symbol")]
t2gsimple <- t2gsimple[-which(duplicated(t2gsimple$ens_gene)),]
t2gsimple$symbol[which(t2gsimple$ens_gene=="ENSRNOG00000003622")]="CYBB"
imexp <- merge(imexp, t2gsimple, by="ens_gene")
imexp$symbol <- toupper(imexp$symbol)
colnames(imexp)[5] <- "gene"
#imexp$gene <- factor(imexp$gene, levels = c("CXCL1", "IFNG", "IL10", "IL13", "IL1B", "IL4", "IL5", "IL6", "TNF"))

comparisons <- list(c("PEpreg", "WTpreg"), c("PEpost", "WTpost"))
imexp$PhenoNames <- factor(imexp$PhenoNames, levels=c("PEpreg", "WTpreg", "PEpost", "WTpost"))

ggboxplot(data=imexp, x="PhenoNames", y="expr", group = "gene", color="PhenoNames", add = "jitter",
          ggtheme = theme_gray()) +
          scale_y_continuous(expand = expansion(mult = .1)) +
          scale_color_manual(values = mcols) +
          stat_compare_means(method = "t.test", comparisons = comparisons,
                     na.rm=T) +
          facet_wrap(~gene, scales="free")

ggsave("../plots/final/hallmark_infl_rna.png", width = 12, height = 8, dpi=150)

# Meso-scale panel

meso <- read.table("../metadata/meso-scale.tsv", sep="\t", header = T)
meso$PhenoNames <- paste(meso$Pheno, meso$Group, sep="")
meso$PhenoNames <- factor(meso$PhenoNames, levels=c("PEpreg", "WTpreg", "PEpost", "WTpost"))
meso$Assay <- factor(meso$Assay, levels = c("KC/GRO", "IFN-gamma", "IL-10", "IL-13", "IL-1b", "IL-4", "IL-5", "IL-6", "TNF-alpha"))

ggboxplot(data=meso, x="PhenoNames", y="Concentration", group = "Assay", color="PhenoNames", add = "jitter",
          ggtheme = theme_gray()) +
          scale_y_continuous(expand = expansion(mult = .1)) +
          scale_color_manual(values = mcols) +
          stat_compare_means(method = "t.test", comparisons = list(c("PEpreg", "WTpreg"), c("PEpost", "WTpost")),
                             na.rm=T) +
          facet_wrap(~Assay, scales="free")
  
ggsave("../plots/final/meso_scale_protein.png", width = 12, height = 8, dpi=150)
  