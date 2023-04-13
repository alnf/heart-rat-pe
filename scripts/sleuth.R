library(sleuth)
library(yaml)
library(dplyr)
library(biomaRt)
library(ggfortify)
library(RColorBrewer)

# Read metadata

metadata <- read.table("../metadata/P2022_metadata.tsv", sep="\t", header = T)
config <- unlist(read_yaml("../config.yaml"))

kal_dirs <- dir(file.path(config["data.kallisto"]))
sample_id <- sapply(strsplit(kal_dirs,"_"), `[`, 1)
sample_id <- sapply(strsplit(sample_id,"-"), `[`, 1)
kal_df <- data.frame(kal_dirs = kal_dirs, SampleID = sample_id)

metadata <- merge(metadata, kal_df, by="SampleID")

# Create transcript to gene mapping
ensembl <- useEnsembl(biomart = "genes", version=109)
ensembl <- useDataset(dataset = "rnorvegicus_gene_ensembl", mart = ensembl)
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Create sleuth objects
summary <- data.frame(sample=metadata$SampleID, condition=metadata$Group,
                 path = paste(config["data.kallisto"], "/", metadata$kal_dirs, sep=""))
summary <- summary %>% filter(grepl("LV", condition))

design <- ~ condition
so <- sleuth_prep(summary, 
                  full_model = design, 
                  target_mapping = t2g, 
                  read_bootstrap_tpm = TRUE,
                  extra_bootstrap_summary = TRUE,
                  transformation_function = function(x) log2(x + 0.5)) 

# Computing model
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

# Check data visually
plot_pca(so, color_by = 'condition')
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
pca <- prcomp(t(log2(sleuth_matrix+0.5)))
autoplot(pca, data = summary, colour = 'condition')

# Getting significant transcripts
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
sleuth_significant <- dplyr::filter(sleuth_significant, var_obs > 10)
head(sleuth_significant)

# Heatmap
my_sample_row <- data.frame(GeneName = sleuth_significant$ext_gene)
rownames(my_sample_row) <- rownames(sleuth_significant)
my_sample_row <- t(my_sample_row)

dev.off()
plot_transcript_heatmap(so, trans="log2", transcripts = sleuth_significant$target_id,
                        labels_row = sleuth_significant$ext_gene,
                        color_low = "#4575B4", color_high = "#D73027", color_mid = "#FFFFFF")


# Gene level statistics
so.g <- sleuth_prep(summary, 
                  full_model = design, 
                  target_mapping = t2g, 
                  read_bootstrap_tpm = TRUE,
                  extra_bootstrap_summary = TRUE,
                  aggregation_column = 'ens_gene',
                  transformation_function = function(x) log2(x + 0.5)) 

so.g <- sleuth_fit(so.g, ~condition, 'full')
so.g <- sleuth_fit(so.g, ~1, 'reduced')
so.g <- sleuth_lrt(so.g, 'reduced', 'full')

sleuth_matrix <- sleuth_to_matrix(so.g, 'obs_norm', 'tpm')
pca <- prcomp(t(log2(sleuth_matrix+0.5)))
autoplot(pca, data = summary, colour = 'condition')

sleuth_table_gene <- sleuth_results(so.g, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)

sleuth_sig_filtered <- sleuth_significant[which(sleuth_significant$ens_gene %in% sleuth_table_gene$target_id[1:10]),]

dev.off()
plot_transcript_heatmap(so, trans="log2", transcripts = sleuth_sig_filtered$target_id,
                        labels_row = sleuth_sig_filtered$ext_gene,
                        color_low = "#4575B4", color_high = "#D73027", color_mid = "#FFFFFF")
