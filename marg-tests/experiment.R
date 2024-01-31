
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

# n_edges <- function(graph) {
#     # ensure symmetric bi-edges
#     graph <- make_symmetric(graph)
#     # now remove them from the lower triangle
#     graph[lower.tri(graph) & graph == 100] <-
#     # return the sum of all non-zero elements: edges
#     return(sum(graph > 0))
# }
n_edges <- function(graph) {
    # represent all edge in upper tri
    graph <- graph + t(graph)
    # how many non-zeros in the upper tri?
    return(length(graph[upper.tri(graph) & graph > 0]))
}

# # Parameter settings
# N_graphs <- 1
# d_nodes <- 4
# max_in_degree <- Inf # default value
# sample_size <- 10000
# fast <- FALSE

run_experiment <- function(
    aridities = c("any", "arid", "maximal-arid", "projection"),
    N_graphs = 10,
    d_nodes = 6,
    max_in_degree = 2,
    sample_size = 10000,
    fast = FALSE,
    verbose = FALSE,
    depth.max = 1, # default
    maxIter = 10, # the default, I think
    faithful.eps = 0, # some default
    equivalent.eps = 1e-10, # some default
    time.max = 10 # default
) {

    n <- sample_size

    # The ground truths already have associated parameters and a covMat
    ground_truths <- list()
    for (i in 1:N_graphs) {
        ground_truths[[i]] <- GenerateAridGroundTruth(
            p = d_nodes,
            max.in.degree = max_in_degree
        )
    }

    empty_graph <- matrix(0, d_nodes, d_nodes)
    colnames(empty_graph) <- rownames(empty_graph) <- colnames(ground_truths[[1]]$mg)

    results <- list()
    results$ground_truths <- ground_truths
    # results$aridity_times <- list()
    true_minabs_causal_effect_mats <- list()
    est_minabs_causal_effect_mats <- list()

    for (aridity in aridities) {
        start_time <- proc.time()
        print(paste("Testing aridity:", aridity))
        results[[aridity]]$true_minabs_causal_effect_mats <- list()
        results[[aridity]]$est_minabs_causal_effect_mats <- list()
        for (i in 1:N_graphs) {
            ground_truth <- ground_truths[[i]]
            # Create data and covariance matrices

            # Can also be handled by the fastGreedySearch function, I think
            data <- GenerateData(n = sample_size, ground_truth$params)
            covMat <- cov(data)

            # Run a single greedy search from the empty graph
            greedy_result <- fastGreedySearch(
                mg.start = empty_graph,
                data=data,
                n=sample_size,
                # maxSteps=Inf, # default
                # direction=3, # stands for "both" forward and backward
                # maxIter=10, # default max RICF iterations
                # edge.penalty=1, # default
                verbose=verbose,
                covMat=covMat,
                # faithful.eps=0, # default, something with the faithfulness test in computeComponentScore, doesn't seems like it's used
                # max.pos=Inf, # default, I think this is the maximum linear index of the graph matrix? not documented
                # dags.only=FALSE, # default,  Don't restrict the search to dags.
                # eps.conv=1e-12, # default, Minimum score difference for greedy search
                aridity=aridity # or "arid", "maximal", "maximal-arid", "projection"
            )

            results[[aridity]][[i]] <- greedy_result

            resulting_graph <- greedy_result$states[[length(greedy_result$states)]]$mg
            resulting_score <- greedy_result$states[[length(greedy_result$states)]]$score

            
            # Find equivalent models of greedy result
            if (fast){
                tmp <- fastFindEquivalentModels(list(mg=resulting_graph), c(), list(),
                                                resulting_score, equivalent.eps, depth.max=depth.max,
                                                n=sample_size, maxIter=maxIter, covMat=covMat,
                                                faithful.eps=faithful.eps)
            } else {
                tmp <- findEquivalentModels(list(mg=resulting_graph), c(), list(),
                                            resulting_score, equivalent.eps, depth.max=depth.max,
                                            n=sample_size, maxIter=maxIter, covMat=covMat,
                                            faithful.eps=faithful.eps)
            }
            scores <- tmp$scores
            equiv.models.greedy <- tmp$res
            print(paste("Number of models in emp. EC of result", i, ":", length(tmp$res)))

            # Find equivalent models of ground truth
            if (fast) {
                tmp <- fastFindEquivalentModels(list(mg=ground_truth$mg), c(), scores, NA, equivalent.eps,
                                                depth.max=depth.max, n=sample_size, maxIter=maxIter,
                                                covMat=covMat, faithful.eps=faithful.eps)
            } else {
                tmp <- findEquivalentModels(list(mg=ground_truth$mg), c(), scores, NA, equivalent.eps,
                                            depth.max=depth.max, time.max=time.max,
                                            n=sample_size, maxIter=maxIter, covMat=covMat, faithful.eps=faithful.eps)
            }
            equiv.models.gt <- tmp$res
            print(paste("Number of models in emp. EC of ground truth:", length(tmp$res)))

            # Compute causal effects of greedy result
            CE.greedy <- lapply(equiv.models.greedy,
                                function(model) getCausalEffects(
                                    fitAncestralGraphCustom(model$mg, covMat, sample_size, maxIter=maxIter)$Bhat))
            CE.minabs.greedy <- do.call("pmin", lapply(CE.greedy, function(mat) abs(mat))) - diag(d_nodes)

            # Compute causal effects of ground truth
            CE.gt <- lapply(equiv.models.gt,
                            function(model) getCausalEffects(
                                fitAncestralGraphCustom(model$mg, covMat, sample_size, maxIter=maxIter)$Bhat))
            CE.minabs.gt <- do.call("pmin", lapply(CE.gt, function(mat) abs(mat))) - diag(d_nodes)

            results[[aridity]]$est_minabs_causal_effect_mats[[i]] <- CE.minabs.greedy
            results[[aridity]]$true_minabs_causal_effect_mats[[i]] <- CE.minabs.gt
            # print("minabs greedy")
            # print(CE.minabs.greedy)
            # print("results minabs greedy")
            # print(results[[aridity]]$est_minabs_causal_effect_mats)
            # return(list(gt=gt, covMat=covMat, res.greedy=res.greedy, i=i, equiv.models.gt=equiv.models.gt,
            #             equiv.models.greedy=equiv.models.greedy, CE.gt=CE.gt, CE.greedy=CE.greedy,
            #             CE.minabs.gt=CE.minabs.gt, CE.minabs.greedy=CE.minabs.greedy, data=data,
            #             n=n))

        }
        end_time <- proc.time()
        elapsed_time <- end_time - start_time
        # results$aridity_times[[aridity]] <- elapsed_time
        results[[aridity]]$aridity_time <- elapsed_time
    }


    # For some reason, extra things are added onto the ends of these.
    # I do not know why.
    # I do know that just removing them fixes it, so that ROC curves can be plotted.
    for(aridity in aridities) {
        results[[aridity]]$true_minabs_causal_effect_mats <- results[[aridity]]$true_minabs_causal_effect_mats[1:N_graphs]
        results[[aridity]]$est_minabs_causal_effect_mats <- results[[aridity]]$est_minabs_causal_effect_mats[1:N_graphs]
    }
    return(results)
}


# thursday night running here

# elapsed: 25039.924 is 6:57:19.924
# start <- proc.time()
# results4in8d100N <- run_experiment(
#     N_graphs = 100,
#     d_nodes = 8,
#     max_in_degree = 4,
#     sample_size = 10000,
#     fast = TRUE,
# )
# end <- proc.time()
# elapsed <- end - start

# saveRDS(results4in8d100N, "./marg-tests/results/exp4in8d100N.RDS")
results4in8d100N <- readRDS("./marg-tests/results/exp4in8d100N.RDS")


results4in8d100N$any$aridity_time / 100
results4in8d100N$arid$aridity_time / 100
results4in8d100N$"maximal-arid"$aridity_time / 100
results4in8d100N$projection$aridity_time / 100

edge_counts4in8d100N <- sapply(results4in8d100N$ground_truths, function(res) n_edges(res$mg))
print(paste("Mean edge counts:", paste(mean(edge_counts4in8d100N))))
print(paste("SD edge counts:", paste(sd(edge_counts4in8d100N))))

results4in8d100N_eval <- list()
options(vsc.dev.args = list(width = 150, height = 150))

pdf("./marg-tests/results/exp4in8d100N_any.pdf", width=7,height=8)
par(mar = c(4,4.5,1,1))
par(cex.axis=1.5, cex.lab=1.5)
results4in8d100N_eval$any <- plotROCCurve(
    results4in8d100N$any$true_minabs_causal_effect_mats,
    results4in8d100N$any$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp4in8d100N_arid.pdf", width =7, height = 8)
par(mar = c(4,4.5,1,1))
par(cex.axis=1.5, cex.lab=1.5)
results4in8d100N_eval$arid <- plotROCCurve(
    results4in8d100N$arid$true_minabs_causal_effect_mats,
    results4in8d100N$arid$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp4in8d100N_maximal-arid.pdf", width =7, height = 8)
par(mar = c(4,4.5,1,1))
par(cex.axis=1.5, cex.lab=1.5)
results4in8d100N_eval$"maximal-arid" <- plotROCCurve(
    results4in8d100N$"maximal-arid"$true_minabs_causal_effect_mats,
    results4in8d100N$"maximal-arid"$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp4in8d100N_projection.pdf", width =7, height = 8)
par(mar = c(4,4.5,1,1))
par(cex.axis=1.5, cex.lab=1.5)
results4in8d100N_eval$projection <- plotROCCurve(
    results4in8d100N$projection$true_minabs_causal_effect_mats,
    results4in8d100N$projection$est_minabs_causal_effect_mats
)
dev.off()

# Get average AUCs and SDs
for (aridity in c("any", "arid", "maximal-arid", "projection")) {
    print(aridity)
    print(paste("Mean", mean(results4in8d100N_eval[[aridity]]$aucs)))
    print(paste("SD:", sd(results4in8d100N_eval[[aridity]]$aucs)))
}
# to here

t.test(
    x = results4in8d100N_eval$any$aucs,
    y = results4in8d100N_eval$projection$aucs,
    paired = TRUE
)
t.test(
    x = results4in8d100N_eval$"maximal-arid"$aucs,
    y = results4in8d100N_eval$projection$aucs,
    # alternative = "greater",
    paired = TRUE
)

# Assuming your data is in a data frame named 'paired_data'
# with columns 'Group1', 'Group2', 'Group3', and 'Group4'


any4 = results4in8d100N_eval$any$aucs
arid4 = results4in8d100N_eval$arid$aucs
max_arid4 = results4in8d100N_eval$"maximal-arid"$aucs
projection4 = results4in8d100N_eval$projection$aucs
qqnorm(any4); qqline(any4, col=2)
qqnorm(arid4); qqline(arid4, col=2)
qqnorm(max_arid4); qqline(max_arid4, col=2)
qqnorm(projection4); qqline(projection4, col=2)

# Create a data frame with paired data
paired_exp4in8d100N <- data.frame(
  any4 = results4in8d100N_eval$any$aucs,
  arid4 = results4in8d100N_eval$arid$aucs,
  max_arid4 = results4in8d100N_eval$"maximal-arid"$aucs,
  projection4 = results4in8d100N_eval$projection$aucs
)

# Add Graph variable for pairing
paired_exp4in8d100N$Graph <- 1:100

# Reshape data to long format
long_exp4in8d100N <- data.frame(
  Graph = rep(paired_exp4in8d100N$Graph, 4),
  Algorithm = rep(colnames(paired_exp4in8d100N[, 1:4]), each = 100),
  Score = c(paired_exp4in8d100N$any4, paired_exp4in8d100N$arid4, paired_exp4in8d100N$max_arid4, paired_exp4in8d100N$projection4)
)
saveRDS(long_exp4in8d100N, "./marg-tests/results/long_exp4in8d100N.RDS")

# Perform repeated measures ANOVA
anova_result <- aov(Score ~ Algorithm + Error(Graph/Algorithm), data = long_data)

# Display the summary
summary(anova_result)


# second dataset

# elapsed: 25776.25 7:09:36.25
# start <- proc.time()
# results3in8d100N <- run_experiment(
#     N_graphs = 100,
#     d_nodes = 8,
#     max_in_degree = 4,
#     sample_size = 10000,
#     fast = TRUE,
# )
# end <- proc.time()
# elapsed3 <- end - start

# saveRDS(results3in8d100N, "./marg-tests/results/exp3in8d100N.RDS")
results3in8d100N <- readRDS("./marg-tests/results/exp3in8d100N.RDS")


results3in8d100N$any$aridity_time / 100
results3in8d100N$arid$aridity_time / 100
results3in8d100N$"maximal-arid"$aridity_time / 100
results3in8d100N$projection$aridity_time / 100


edge_counts3in8d100N <- sapply(results3in8d100N$ground_truths, function(res) n_edges(res$mg))
print(paste("Mean edge counts:", paste(mean(edge_counts3in8d100N))))
print(paste("SD edge counts:", paste(sd(edge_counts3in8d100N))))

options(vsc.dev.args = list(width = 250, height = 300))
results3in8d100N_eval <- list()
pdf("./marg-tests/results/exp3in8d100N_any.pdf")
results3in8d100N_eval$any <- plotROCCurve(
    results3in8d100N$any$true_minabs_causal_effect_mats,
    results3in8d100N$any$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp3in8d100N_arid.pdf")
results3in8d100N_eval$arid <- plotROCCurve(
    results3in8d100N$arid$true_minabs_causal_effect_mats,
    results3in8d100N$arid$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp3in8d100N_maximal-arid.pdf")
results3in8d100N_eval$"maximal-arid" <- plotROCCurve(
    results3in8d100N$"maximal-arid"$true_minabs_causal_effect_mats,
    results3in8d100N$"maximal-arid"$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp3in8d100N_projection.pdf")
results3in8d100N_eval$projection <- plotROCCurve(
    results3in8d100N$projection$true_minabs_causal_effect_mats,
    results3in8d100N$projection$est_minabs_causal_effect_mats
)
dev.off()


# Get average AUCs and SDs
for (aridity in c("any", "arid", "maximal-arid", "projection")) {
    print(aridity)
    print(paste("Mean", mean(results3in8d100N_eval[[aridity]]$aucs)))
    print(paste("SD:", sd(results3in8d100N_eval[[aridity]]$aucs)))
}

paired_exp3in8d100N <- data.frame(
  any3 = results3in8d100N_eval$any$aucs,
  arid3 = results3in8d100N_eval$arid$aucs,
  max_arid3 = results3in8d100N_eval$"maximal-arid"$aucs,
  projection3 = results3in8d100N_eval$projection$aucs
)

# Add Graph variable for pairing
paired_exp3in8d100N$Graph <- 1:100

# Reshape data to long format
long_exp3in8d100N <- data.frame(
  Graph = rep(paired_exp3in8d100N$Graph, 4),
  Algorithm = rep(colnames(paired_exp3in8d100N[, 1:4]), each = 100),
  Score = c(paired_exp3in8d100N$any3, paired_exp3in8d100N$arid3, paired_exp3in8d100N$max_arid3, paired_exp3in8d100N$projection3)
)
saveRDS(long_exp3in8d100N, "./marg-tests/results/long_exp3in8d100N.RDS")

# to here

t.test(
    x = results3in8d100N_eval$any$aucs,
    y = results3in8d100N_eval$projection$aucs,
    paired = TRUE
)
t.test(
    x = results3in8d100N_eval$"maximal-arid"$aucs,
    y = results3in8d100N_eval$projection$aucs,
    alternative = "less",
    paired = TRUE
)





# Yet another one


# start <- proc.time()
# results3in6d200N <- run_experiment(
#     max_in_degree = 3,
#     d_nodes = 6,
#     N=200,
#     fast=TRUE,
# )
# end <-proc.time()
# elapsed2in6d200N <- end - start
# elapsed2in6d200N
# saveRDS(results3in6d200N, "./marg-tests/results/results3in6d200N.RDS")
results3in6d200N <- readRDS("./marg-tests/results/results3in6d200N.RDS")

results3in6d200N_eval <- list()
for (aridity in c("any", "arid", "maximal-arid", "projection")) {
    pdf(paste0("./marg-tests/results/exp3in6d100N", aridity, ".pdf"))
    results3in6d200N_eval[[aridity]] <- plotROCCurve(
        results3in6d200N[[aridity]]$true_minabs_causal_effect_mats,
        results3in6d200N[[aridity]]$est_minabs_causal_effect_mats
    )
    dev.off()
}


results3in6d200N$any$aridity_time / 100
results3in6d200N$arid$aridity_time / 100
results3in6d200N$"maximal-arid"$aridity_time / 100
results3in6d200N$projection$aridity_time / 100


edge_counts3in6d200N <- sapply(results3in6d200N$ground_truths, function(res) n_edges(res$mg))
print(paste("Mean edge counts:", mean(edge_counts3in6d200N)))
print(paste("SD edge counts:", sd(edge_counts3in6d200N)))

# Get average AUCs and SDs
for (aridity in c("any", "arid", "maximal-arid", "projection")) {
    print(aridity)
    print(paste("Mean", mean(results3in6d200N_eval[[aridity]]$aucs)))
    print(paste("SD:", sd(results3in6d200N_eval[[aridity]]$aucs)))
}

paired_exp3in6d200N <- data.frame(
  any = results3in6d200N_eval$any$aucs,
  arid = results3in6d200N_eval$arid$aucs,
  max_arid = results3in6d200N_eval$"maximal-arid"$aucs,
  projection = results3in6d200N_eval$projection$aucs
)

# Add Graph variable for pairing
paired_exp3in6d200N$Graph <- 1:200

# Reshape data to long format
long_exp3in6d200N <- data.frame(
  Graph = rep(paired_exp3in6d200N$Graph, 4),
  Algorithm = rep(colnames(paired_exp3in6d200N[, 1:4]), each = 200),
  Score = c(paired_exp3in6d200N$any, paired_exp3in6d200N$arid, paired_exp3in6d200N$max_arid, paired_exp3in6d200N$projection)
)
saveRDS(long_exp3in6d200N, "./marg-tests/results/long_exp3in6d200N.RDS")




# Never properly run:
# start <- proc.time()
# results4in10d100N <- run_experiment(
#     N_graphs = 10,
#     d_nodes = 5,
#     max_in_degree = 3,
#     sample_size = 10000,
#     fast = TRUE,
#     aridities = c("projection")
# )
# end <- proc.time()
# elapsed <- end - start

# saveRDS(results4in10d100N, "./marg-tests/results/exp4in10d100N.RDS")

# options(vsc.dev.args = list(width = 250, height = 300))

# pdf("./marg-tests/results/exp4in10d100N_any.pdf")
# plotROCCurve(
#     results4in10d100N$any$true_minabs_causal_effect_mats,
#     results4in10d100N$any$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp4in10d100N_arid.pdf")
# plotROCCurve(
#     results4in10d100N$arid$true_minabs_causal_effect_mats,
#     results4in10d100N$arid$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp4in10d100N_maximal-arid.pdf")
# plotROCCurve(
#     results4in10d100N$"maximal-arid"$true_minabs_causal_effect_mats,
#     results4in10d100N$"maximal-arid"$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp4in10d100N_projection.pdf")
# plotROCCurve(
#     results4in10d100N$projection$true_minabs_causal_effect_mats,
#     results4in10d100N$projection$est_minabs_causal_effect_mats
# )
# dev.off()

# for(aridity in c("any", "arid", "maximal-arid", "projection")) {
#     pdf(paste0("./marg-tests/results/exp4in10d100N_", aridity, ".pdf"))
#     plotROCCurve(
#         results4in10d100N[[aridity]]$true_minabs_causal_effect_mats,
#         results4in10d100N[[aridity]]$est_minabs_causal_effect_mats
#     )
#     dev.off()
# }



# unused
# ### Tests of N=100 (d=6 max in =2)
# # not reported

# # start <- proc.time()
# # results100 <- run_experiment(N=100, fast=TRUE)
# # end <-proc.time()
# # elapsed2in6d100N <- end - start
# # elapsed2in6d100N
# # saveRDS(results100, "./marg-tests/results/results100_v2.RDS")
# results100 <- readRDS("./marg-tests/results/results100_v2.RDS")


# options(vsc.dev.args = list(width = 250, height = 300))
# pdf("./marg-tests/results/exp100any.pdf")
# plotROCCurve(
#     results100[["any"]]$true_minabs_causal_effect_mats,
#     results100[["any"]]$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp100arid.pdf")
# plotROCCurve(
#     results100[["arid"]]$true_minabs_causal_effect_mats,
#     results100[["arid"]]$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp100maximal-arid.pdf")
# plotROCCurve(
#     results100[["maximal-arid"]]$true_minabs_causal_effect_mats,
#     results100[["maximal-arid"]]$est_minabs_causal_effect_mats
# )
# dev.off()
# pdf("./marg-tests/results/exp100maximal-arid.pdf")
# plotROCCurve(
#     results100[["projection"]]$true_minabs_causal_effect_mats,
#     results100[["projection"]]$est_minabs_causal_effect_mats
# )
# dev.off()

# start <- proc.time()
# results100projection <- run_experiment(N=100, aridities="projection", fast=TRUE)
# end <-proc.time()
# elapsed <- end - start
# elapsed
# saveRDS(results100projection, "./marg-tests/results/results100projection.RDS")

# pdf("./marg-tests/results/exp100projection.pdf")
# plotROCCurve(
#     results100projection[["projection"]]$true_minabs_causal_effect_mats,
#     results100projection[["projection"]]$est_minabs_causal_effect_mats
# )
# dev.off()

# # Run a single greedy search from the empty graph
# greedy_result <- fastGreedySearch(
#     mg.start = empty_graph,
#     # data=NULL,
#     n=sample_size,
#     # maxSteps=Inf, # default
#     # direction=3, # stands for "both" forward and backward
#     # maxIter=10, # default max RICF iterations
#     # edge.penalty=1, # default
#     # verbose=TRUE, # default
#     # covMat=NULL,
#     # faithful.eps=0, # default, something with the faithfulness test in computeComponentScore, doesn't seems like it's used
#     # max.pos=Inf, # default, I think this is the maximum linear index of the graph matrix? not documented
#     # dags.only=FALSE, # default,  Don't restrict the search to dags.
#     # eps.conv=1e-12, # default, Minimum score difference for greedy search
#     aridity=aridity # or "arid", "maximal", "maximal-arid", "projection"
# )

for (state in greedy_result$states) {
    graph <- state$mg
    if (!all(graph == 0)) {
        print(graph)
        plotBAP(graph)
        Sys.sleep(3)
    } else {
        print(graph)
        Sys.sleep(3)
    }
}
graphs <- lapply(greedy_result$states, function(s) s$mg)
