library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(rbioapi)
library(STRINGdb)
library(cowplot)
library(patchwork)

# String IDs
modname = "m2"
rat_ids <- rat[which(rat$modules == modname),]
rat_ids$symbol[which(rat_ids$symbol=="")]=rat_ids$msymbol[which(rat_ids$symbol=="")]
rat_ids <- rat_ids[-which(duplicated(rat_ids$ens_gene)),]
str_ids <- rba_string_map_ids(ids = rat_ids$symbol, species = 10116)
fun_ids <- ora[which(ora$Module==toupper(modname)),]$geneID

# 1. Read the edge list
#edges <- read.table("../metadata/string_m2.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#nodes <- read.table("../metadata/string_m2_nodes.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

edges <- rba_string_interactions_network(str_ids$stringId, species = 10116)

# 2. Create an igraph graph object
g <- igraph::graph_from_data_frame(edges[3:ncol(edges)], directed=FALSE)

# 3. Optionally convert igraph to a tidygraph object for easy plotting with ggraph
g_tidy <- as_tbl_graph(g)

# 4. Add zero degree nodes
num_nodes <- g_tidy %>% 
  as_tibble("nodes") %>% 
  nrow()

#znodes <- tibble(name = nodes$node[which(nodes$node_degree==0)])

#g_tidy <- g_tidy %>%
#  bind_nodes(znodes)

# 5. Visualize with ggraph
fixed_layout <- create_layout(g_tidy, layout = "fr")

#################
plot_list <- c()
for (i in 1:length(fun_ids)) {
  ll <- list(ens_gene = fun_ids[[i]])
  ids <- merge(ll, rat_ids, by="ens_gene")
  print(ids$symbol)
  pp <- ggraph(fixed_layout) +
    geom_edge_link(alpha = 0.5) +
    geom_node_point(aes(color = name %in% ids$symbol), size = 3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    theme_graph() +
    theme(legend.position="none", plot.title = element_text(size=11)) +
    ggtitle(names(fun_ids[i]))
  pp
  plot_list <- c(plot_list, list(pp))
}


plot_list <- lapply(seq_along(fun_ids), function(i) {
  # 1) Subset/merge your IDs:
  df_ids <- data.frame(ens_gene = fun_ids[[i]])
  ids <- merge(df_ids, rat_ids, by = "ens_gene")
  
  # 2) Build the plot in this local function environment:
  ggraph(fixed_layout) +
    geom_edge_link(alpha = 0.5) +
    geom_node_point(aes(color = name %in% ids$symbol), size = 3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    theme_graph() +
    theme(legend.position="none", plot.title = element_text(size=11)) +
    ggtitle(names(fun_ids[i]))
})

wp <- wrap_plots(plot_list, ncol=3)
ggsave(filename = "../plots/final/string.png", plot = wp, dpi=100, width=18, height=27)
