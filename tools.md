# Tools

Main software and versions used in the project:

- Reference genome: mRatBN7.2 (GCA_015227675.2)
- Quality control: fastp 0.22.0 (default parameters)
- Quantification: kallisto 0.48.0 (default parameters)
- Differential expression: DESeq2 1.38.3
- Visualisation: ggplot2 3.4.4, pheatmap 1.0.12, VennDiagram 1.7.3
- Enrichment analysis: ssGSEA2 1.0.0, cmapR 1.10.0, ShinyGO 0.80
- R version 4.2.3
- BioConductor 3.16
- BioMart 109

## Environments


Conda environment export:

```
conda env export > ratseq.yml
```

R environment session information:

```
sink("sessionInfo.txt")
sessionInfo()
sink()
```

You can find both files in the root of the repository.