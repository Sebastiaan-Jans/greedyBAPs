
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/api.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/greedy.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/icfmagCustom.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/plotting.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/simulation.R", encoding = "UTF-8")

source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/margs.R", encoding = "UTF-8")
library(foreach) # for %dopar% (paralellism)
library(ggm)

# What I *should* do:
# library(greedyBAPs)

set.seed(1)
# How many d-noded uniformly sampled BAPs are (maximal) arid?

N_graphs <- 100
d_nodes <- 10
max_in_degree <- Inf

arid_count <- 0
for (i in 1:N_graphs) {
    arid <- is_arid(GenerateMG4(p = d_nodes, max.in.degree = max_in_degree)[[1]])
    print(arid)
    arid_count <- arid_count + arid
}
arid_count / N_graphs

aridity_heatmap <- function(
    N_graphs = 10,
    model_sizes = c(2, 3, 4, 5, 6, 7, 8, 9, 10),
    max_in_degrees = c(1, 2, 3, 4, 5, Inf),
    aridity_func = is_arid
) {
    set.seed(1)
    heatmap <- matrix(0, nrow = length(max_in_degrees), ncol = length(model_sizes))
    colnames(heatmap) <- model_sizes
    rownames(heatmap) <- max_in_degrees
    for (j in 1:length(model_sizes)) {
        d_nodes <- model_sizes[j]
        for (i in 1:length(max_in_degrees)) {
            max_in_degree <- max_in_degrees[i]
            count <- 0
            for (k in 1:N_graphs) {
                arid <- aridity_func(GenerateMG4(p = d_nodes, max.in.degree = max_in_degree)[[1]])
                count <- count + arid
            }
            heatmap[i, j] <- count/N_graphs
            print(paste0("With ", d_nodes, " nodes and max in-degree ", max_in_degree, ", ", 100*count/N_graphs, "%"))
        }
    }
    return(heatmap)
}

# Load ggplot2 package
library(ggplot2)
library(reshape)

# arid_start <- proc.time()
# arid_heatmap <- 100*aridity_heatmap(N=1000, aridity_func = is_arid)
# arid_end <- proc.time()
# arid_elapsed <- arid_end - arid_start
# arid_elapsed

# saveRDS(arid_heatmap, file="./marg-tests/results/arid_heatmap_1000.RDS")
arid_heatmap_big <- readRDS("./marg-tests/results/arid_heatmap_1000.RDS")
arid_heatmap <- arid_heatmap_big[, 2:ncol(arid_heatmap_big)]

# Reshape data
arid_melt <- melt(arid_heatmap)
colnames(arid_melt) <- c("max_in_degree", "model_size", "value")
arid_melt$max_in_degree <- as.character(arid_melt$max_in_degree)

options(vsc.dev.args = list(width = 400, height = 200))
# Create heatmap
arid_plot <- ggplot(arid_melt, aes(x = model_size, y = max_in_degree)) +
    geom_tile(aes(fill = value)) +
    labs(
        x = "Model size (nodes)",
        y = "Max in-degree",
        # title = "Percentage of BAPs that are maximal-arid",
        fill = "% arid"
    ) +
    scale_x_discrete(limits=c(3,4,5,6,7,8,9,10)) +
    geom_text(aes(x= model_size, y = max_in_degree, label=round(value)), color="white") +
    scale_fill_continuous(limits = c(0, 100))
arid_plot

ggsave("./marg-tests/results/arid_heatmap.pdf", plot = arid_plot, device = "pdf", width = 5, height = 2.5)

# maximal_start <- proc.time()
# maximal_arid_heatmap <- 100 * aridity_heatmap(N=1000, aridity_func = is_maximal_arid)
# maximal_end <- proc.time()
# maximal_elapsed <- maximal_end - maximal_start
# maximal_elapsed

# saveRDS(maximal_arid_heatmap, "./marg-tests/results/maximal_arid_heatmap_1000.RDS")
maximal_arid_heatmap_big <- readRDS("./marg-tests/results/maximal_arid_heatmap_1000.RDS")
maximal_arid_heatmap <- maximal_arid_heatmap_big[, 2:ncol(maximal_arid_heatmap_big)]

# Reshape data
max_melt <- melt(maximal_arid_heatmap)
colnames(max_melt) <- c("max_in_degree", "model_size", "value")
max_melt$max_in_degree <- as.character(max_melt$max_in_degree)

options(vsc.dev.args = list(width = 400, height = 200))
# Create heatmap
maximal_arid_plot <- ggplot(max_melt, aes(x = model_size, y = max_in_degree)) +
    geom_tile(aes(fill = value)) +
    labs(
        x = "Model size (nodes)",
        y = "Max in-degree",
        # title = "Percentage of BAPs that are maximal-arid",
        fill = "% MArG"
    ) +
    scale_x_discrete(limits=c(3,4,5,6,7,8,9,10)) +
    geom_text(aes(x= model_size, y = max_in_degree, label=round(value)), color="white") +
    scale_fill_continuous(limits = c(0, 100))
maximal_arid_plot

ggsave("./marg-tests/results/maximal_arid_heatmap.pdf", plot = maximal_arid_plot, device = "pdf", width = 5, height = 2.5)




# ================== Some extra unused stuff ====================
# max_arid_small <- aridity_heatmap(
#     N=1000,
#     model_sizes = c(2, 3),
#     max_in_degrees = c(2, 3, 4, 5, Inf),
#     aridity_func = is_maximal_arid
# )
# max_arid_small2 <- aridity_heatmap(
#     N=100,
#     model_sizes = c(2, 3, 4, 5, 6),
#     max_in_degrees = c(1, 2),
#     aridity_func = is_maximal_arid
# )

# # Do 3-node graphs exists that are non-maximal?
# get_non_max_arid_graphs <- function(N=100, nodes=5, in_degree=2) {
#     graphs <- list()
#     ind <- 1
#     for (i in 1:N) {
#         graph <- GenerateMG4(p=nodes, max.in.degree = in_degree)[[1]]
#         if (!is_maximal_arid(graph)) {
#             graphs[[ind]] <- graph
#             ind <- ind + 1
#             # graphs <- c(graphs, graph)
#         }
#     }
#     return(graphs)
# }