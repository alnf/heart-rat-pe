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
             #"ENSRNOG00000016308", #IL-10ra
             #"ENSRNOG00000028638", #IL-10rb
             "ENSRNOG00000007652", #IL-13
             "ENSRNOG00000002802", #CXCL-1
             "ENSRNOG00000070745"  #TNF
             )

impanel <- c("ENSRNOG00000004647", #IL-10
             "ENSRNOG00000016308", #IL-10ra
             "ENSRNOG00000028638"  #IL-10rb
            )

impanel <- c("ENSRNOG00000010208", #TIMP1
             "ENSRNOG00000004273", #IFITM1
             "ENSRNOG00000003622", #CYBB
             "ENSRNOG00000007159", #CCL2
             "ENSRNOG00000017819", #CD14
             "ENSRNOG00000000239", #CCL7
             "ENSRNOG00000012779", #MSR1
             "ENSRNOG00000014361", #EDN1
             "ENSRNOG00000056219", #OLR1
             "ENSRNOG00000004500", #MYC
             "ENSRNOG00000046254", #ADGRE1
             "ENSRNOG00000004380", #IL12B
             "ENSRNOG00000005479", #SLC1A2
             "ENSRNOG00000016278", #CCL17
             "ENSRNOG00000010665", #CCR7
             "ENSRNOG00000065741"  #IL7R
            )

impanel = c("ENSRNOG00000001959",
            "ENSRNOG00000017414",
            "ENSRNOG00000001963",
            "ENSRNOG00000022839",
            "ENSRNOG00000069835",
            "ENSRNOG00000028768",
            "ENSRNOG00000056947",
            "ENSRNOG00000001187",
            "ENSRNOG00000037198",
            "ENSRNOG00000049282",
            "ENSRNOG00000001963"
           )

impanel = infs
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
imexp$expr <- imexp$expr + 0.001

imexp$PhenoNames <- factor(imexp$PhenoNames)
comparisons <- list(c("WTpreg", "WTpost"), c("PEpreg", "PEpost"), c("WTpreg", "PEpreg"), c("WTpost", "PEpost"))

ggboxplot(data=imexp, x="PhenoNames", y="expr", color="PhenoNames", add = "jitter") +
          stat_compare_means(method = "t.test", comparisons = comparisons) +
          scale_color_manual(values=annoCol$Group) +
          facet_wrap(~symbol, ncol=6)

ggsave("../plots/pe_preg/immuno_panel.png", width = 5, height = 30, dpi = 80)


tpm_abd[which(rownames(tpm_abd)=="ENSRNOG00000010208"),]
lcounts[which(rownames(lcounts)=="ENSRNOG00000010208"),]


# Heatmap

m <- lcounts[impanel, which(colnames(lcounts) %in% metadata.left$SampleNumber)]
m$ens_gene <- rownames(m)
m <- merge(m, t2gh, by="ens_gene")
m <- m[!duplicated(m$hsymbol),]
rownames(m) <- m$hsymbol
m <- m[-which(m$hsymbol=="IFITM2"),]
m <- m[-which(m$hsymbol=="IFITM3"),]
m <- as.matrix(m[,2:(ncol(m)-1)])

m <- m[,match(metadata.left$SampleNumber, colnames(m))]
ha = HeatmapAnnotation(Group = metadata.left[,c("PhenoNames")], col=annoCol, show_legend = F, show_annotation_name = F)

colors = hcl.colors(30, palette = "Purple-Green")
ht = Heatmap(t(scale(t(m))), show_row_names = T, show_row_dend = F, show_column_names = T, cluster_columns = F,
             top_annotation = ha, name = "expr", column_split = cspl)
png("../plots/pe_preg/heatmap_main_genes.png", width = 10.5, height = 19, units="in", res=80)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()


# Meso-scale panel

meso <- read.table("../metadata/meso-scale.tsv", sep="\t", header = T)
meso$PhenoNames <- paste(meso$Pheno, meso$Group, sep="")
meso$PhenoNames <- factor(meso$PhenoNames, levels=c("PEpreg", "WTpreg", "PEpost", "WTpost"))

ggboxplot(data=meso, x="PhenoNames", y="Concentration", group = "Assay", color="PhenoNames", add = "jitter",
          ggtheme = theme_gray()) +
          scale_y_continuous(expand = expansion(mult = .1)) +
          scale_color_manual(values = mcols) +
          stat_compare_means(method = "t.test", comparisons = list(c("PEpreg", "WTpreg"), c("PEpost", "WTpost")),
                             na.rm=T) +
          facet_wrap(~Assay, scales="free")
  
  
  