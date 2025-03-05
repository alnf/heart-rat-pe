library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(rbioapi)
library(STRINGdb)
library(cowplot)
library(patchwork)

# String IDs
#rat_ids <- rat
rat_ids <- rat[which(rat$modules %in% c("m1", "m2")),]

rat_ids$symbol[which(rat_ids$symbol=="")]=rat_ids$msymbol[which(rat_ids$symbol=="")]
rat_ids <- rat_ids[-which(duplicated(rat_ids$ens_gene)),]
str_ids <- rba_string_map_ids(ids = rat_ids$symbol, species = 10116)
#rat_ids <- rat_ids[-which(duplicated(rat_ids$symbol)),]
rat_ids <- merge(rat_ids, str_ids, by.x="symbol", by.y="queryItem")
rat_ids <- merge(rat_ids, PEpp_SDpp[,c("ens_gene", "log2FoldChange")], by="ens_gene")
fun_ids <- ora$geneID

# 1. Read the edge list
#edges <- read.table("../metadata/string_m2.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#nodes <- read.table("../metadata/string_m2_nodes.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

edges <- rba_string_interactions_network(str_ids$stringId, species = 10116, required_score=700)

# 2. Create an igraph graph object
g <- igraph::graph_from_data_frame(edges[3:ncol(edges)], directed=FALSE)

# 3. Optionally convert igraph to a tidygraph object for easy plotting with ggraph
g_tidy <- as_tbl_graph(g)

# 4. Add zero degree nodes
num_nodes <- g_tidy %>% 
  as_tibble("nodes") %>% 
  nrow()
num_nodes

#znodes <- tibble(name = nodes$node[which(nodes$node_degree==0)])

#g_tidy <- g_tidy %>%
#  bind_nodes(znodes)

# 5. Visualize with ggraph
g_tidy_merged <- g_tidy %>%
  mutate(name = as.character(name)) %>%
  left_join(rat_ids, by = c("name" = "preferredName"))

fixed_layout <- create_layout(g_tidy_merged, layout = "fr")
nrow(fixed_layout)

nodcols <- mcols[9:14]
names(nodcols) <- c("m1", "m2", "m3", "m4", "m5", "NA")

# Plot all modules
col_fun = colorRamp2(c(min(rat_ids$log2FoldChange), 0, max(rat_ids$log2FoldChange)), c("blue", "white", "red"))

pp <- ggraph(fixed_layout) +
  geom_edge_link(alpha = 0.5, edge_colour="black") +
  geom_node_point(aes(color = modules, fill=sign(log2FoldChange), size=abs(log2FoldChange)), shape = 21, stroke = 1.5) +
  scale_color_manual(values = nodcols) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_graph()
pp

# Another layout

g_tidy_main_comp <- g_tidy_merged %>%
  mutate(component = group_components()) %>%
  filter(component == which.max(table(component)))

g_tidy_main_comp <- g_tidy_main_comp %>%
  mutate(betweenness = centrality_betweenness()) %>%
  mutate(deg = centrality_degree()) %>%
  mutate(clo = centrality_closeness())

layout_centr <- create_layout(g_tidy_main_comp, layout = "centrality", cent = betweenness)

pp <- ggraph(layout_centr) +
  geom_edge_link(alpha = 0.5, edge_colour="black") +
  geom_node_point(aes(color = modules, size = betweenness)) +
  scale_color_manual(values = nodcols) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_graph()
pp

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
