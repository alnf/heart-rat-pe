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


## Plotting heatmaps

### preg vs post

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

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, csplit=T, dlists, cols, acols, clegend=F, clcol=F, rsplit=rspl$modules, modcols)
png("../plots/final/final_heatmap_genes_pe_reduced_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

### preg, post, preg_effect

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

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, csplit=T, dlists, cols, acols, clegend=F, clcol=F, rsplit=rspl$modules, modcols)
png("../plots/final/final_heatmap_genes_pe_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()






### preg vs post
names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/final_venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost")
makeVenn(2, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)]),
         names,
         paste("DEGs in pregnancy and postpartum,\n lFC=0.58, ngenes=", nrow(genes_pe_reduced), " (genes_pe_reduced list)", sep=""),
         fname, brewer.pal(3, "Pastel2")[1:2], width=600, height=550)





pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp)
cols <- list(preg="pink", post="lightblue")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_reduced.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=T, dlists, cols, annoCol, clegend=F)
png("../plots/heatmap_genes_pe_reduced_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_reduced_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

### preg, post, preg_effect
names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, ".png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost", "PEpreg_PEpost", "WTpreg_WTpost")
makeVenn(4, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)],
                 PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene),
         names,
         paste("DEGs in pregnancy and postpartum and pregnancy effect,\n lFC=0.58, ngenes=", nrow(genes_pe), " (genes_pe list)", sep=""),
         fname, brewer.pal(4, "Pastel2"), dist=0.15, width=800, height=650)


pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp, pe_pref=PEd21_PEpp, wt_pref=SDd21_SDpp)
cols <- list(preg="pink", post="lightblue", pe_pref="orange", wt_pref="maroon")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=T, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe, csplit=T, dlists, cols, annoCol, clegend=T, clcol=F)
png("../plots/heatmap_genes_pe_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

### norm pregnancy
names = c("SDd21_SDpp", "SDd21_np", "SDpp_np")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, ".png", sep="")
names = c("WTpreg_WTpost", "WTpreg_WTnp", "WTpost_WTnp")
makeVenn(3, list(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene), names,
         paste("DEGs in normal pregnancy,\n lFC=0.58, ngenes=", nrow(genes_norm), " (genes_norm list)", sep=""),
         fname, brewer.pal(3, "Pastel2"))



## Plotting venns for lFC=1

### preg vs post

names = c("PEd21_SDd21" , "PEpp_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost")
makeVenn(2, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)]),
         names,
         paste("DEGs in pregnancy and postpartum,\n lFC=1, ngenes=", nrow(genes_pe_reduced_lfc1), " (genes_pe_reduced_lfc1 list)", sep=""),
         fname, brewer.pal(3, "Pastel2")[1:2], width=600, height=550)

pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp)
cols <- list(preg="pink", post="lightblue")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_reduced_lfc1.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=T, dlists, cols, annoCol, clegend=F)
png("../plots/heatmap_genes_pe_reduced_lfc1_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_reduced_lfc1_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_reduced_lfc1, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_reduced_lfc1_ordered_fl.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(1, 1, 1, 2), "mm"))
dev.off()

library(InteractiveComplexHeatmap)
htShiny(ht)


#### Modules
m1 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M1"),]$genes)
m2 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M2"),]$genes)
m3 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M3"),]$genes)

IFN <- data.frame(ens_gene = infs)

dlists <- list(preg=PEd21_SDd21, pregsig=d21, post=PEpp_SDpp, postsig=pp, M1=m1, M2=m2, M3=m3, IFN = IFN, )
cols <- list(preg="pink", pregsig="#585050", post="lightblue", postsig="darkgray",
             pe_pref="orange", wt_pref="maroon",
             M1="#820D3F", M2="#E64A00", M3="#3B3EDE", M4="#871C9A", IFN = "darkgreen", Not.Correlated="gray")

rspl <- data.frame(modules = cem@module$modules)
rspl$ens_gene <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M5", "M3", "M4", "NC"))

glist <- genes_pe_reduced_lfc1
glist$ens_gene <- genes_pe_reduced_lfc1$ens_gene[order(match(genes_pe_reduced_lfc1$ens_gene, rspl$ens_gene))]

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F, rsplit=rspl$modules)
png("../plots/heatmap_genes_pe_reduced_lfc1_ordered_modules_ifn_fl.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

### preg, post, preg_effect
names = c("PEd21_SDd21" , "PEpp_SDpp", "PEd21_PEpp", "SDd21_SDpp")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("PEpreg_WTpreg", "PEpost_WTpost", "PEpreg_PEpost", "WTpreg_WTpost")
makeVenn(4, list(PEd21_SDd21$ens_gene[-which(PEd21_SDd21$ens_gene==agt$ens_gene)], PEpp_SDpp$ens_gene[-which(PEpp_SDpp$ens_gene==agt$ens_gene)],
                 PEd21_PEpp$ens_gene, SDd21_SDpp$ens_gene),
         names,
         paste("DEGs in pregnancy and postpartum and pregnancy effect,\n lFC=1, ngenes=", nrow(genes_pe_lfc1), " (genes_pe_lfc1 list)", sep=""),
         fname, brewer.pal(4, "Pastel2"), dist=0.15, width=800, height=650)

pheno <- c("PEd21", "SDd21", "PEpp", "SDpp")
region <- "LV"
order = NULL
dlists <- list(preg=PEd21_SDd21, post=PEpp_SDpp, pe_pref=PEd21_PEpp, wt_pref=SDd21_SDpp)
cols <- list(preg="pink", post="lightblue", pe_pref="orange", wt_pref="maroon")

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=F, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_lfc1.png", width = 10.5, height = 9.1, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=T, dlists, cols, annoCol, clegend=T)
png("../plots/heatmap_genes_pe_lfc1_csplit.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

order <- c("WTpreg", "WTpost", "PEpreg", "PEpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=T, dlists, cols, annoCol, clegend=T, clcol=F)
png("../plots/heatmap_genes_pe_lfc1_ordered.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, genes_pe_lfc1, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F)
png("../plots/heatmap_genes_pe_lfc1_ordered_fl.png", width = 10.5, height = 9.1, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

#### Modules
m1 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M1"),]$genes)
m2 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M2"),]$genes)
m3 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M3"),]$genes)
m4 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M4"),]$genes)
m5 <- data.frame(ens_gene=cem@module[which(cem@module$modules=="M5"),]$genes)

IFN <- data.frame(ens_gene = infs)
E2F <- data.frame(ens_gene = e2f)
MIS <- data.frame(ens_gene = mis)
G2M <- data.frame(ens_gene = g2m)
dlists <- list(preg=PEd21_SDd21,  pregsig=d21, post=PEpp_SDpp, postsig=pp,
               pe_pref=PEd21_PEpp, wt_pref=SDd21_SDpp,
               M1=m1, M2=m2, M3=m3, M4=m4, M5=m5,
               IFN=IFN, E2F=E2F, MIS=MIS, G2M=G2M)
cols <- list(preg="pink", pregsig="#585050", post="lightblue", postsig="darkgray", 
             pe_pref="orange", wt_pref="maroon",
             M1="#820D3F", M2="#E64A00", M3="#3B3EDE", M4="#871C9A", M5="#14C7BA",
             IFN="darkgreen", E2F="green", MIS="maroon", G2M="darkblue", Not.Correlated="gray")

rspl <- data.frame(modules = cem@module$modules)
rspl$ens_gene <- rownames(m)
rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M5", "M3", "M4", "NC"))

glist <- genes_pe_lfc1
glist$ens_gene <- genes_pe_lfc1$ens_gene[order(match(genes_pe_lfc1$ens_gene, rspl$ens_gene))]

order <- c("PEpreg", "WTpreg", "PEpost", "WTpost")
ht <- compHeatmap(pheno, region, order, metadata.h, lcounts, glist, csplit=T, dlists, cols, annoCol, clegend=F, clcol=F, rsplit=rspl$modules)
png("../plots/heatmap_genes_pe_lfc1_ordered_modules_ifn_fl.png", width = 10.5, height = 9, units="in", res=100)
draw(ht, padding = unit(c(2, 2, 2, 3), "mm"))
dev.off()

### norm pregnancy
names = c("SDd21_SDpp", "SDd21_np", "SDpp_np")
fname = paste(names, collapse="_")
fname = paste("../plots/venn_", fname, "_lfc1.png", sep="")
names = c("WTpreg_WTpost", "WTpreg_WTnp", "WTpost_WTnp")
makeVenn(3, list(SDd21_SDpp$ens_gene, SDd21_np$ens_gene, SDpp_np$ens_gene), names,
         paste("DEGs in normal pregnancy,\n lFC=1, ngenes=", nrow(genes_norm_lfc1), " (genes_norm_lfc1 list)", sep=""),
         fname, brewer.pal(3, "Pastel2"), width=800, height=750)
