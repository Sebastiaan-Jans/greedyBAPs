
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
    true_minabs_causal_effect_mats <- list()
    est_minabs_causal_effect_mats <- list()

    for (aridity in aridities) {
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
            print(paste("Number of models in emp. EC of result:", length(tmp$res)))

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
    }

    # plotROCCurve(
    #     true_minabs_causal_effect_mats,
    #     est_minabs_causal_effect_mats
    # )
    return(results)
}

start <- proc.time()
results100 <- run_experiment(N=100, fast=TRUE)
end <-proc.time()
elapsed <- end - start
elapsed
saveRDS(results100, "./marg-tests/results/results100.RDS")
results100 <- readRDS("./marg-tests/results/results100.RDS")

plotROCCurve(
    results[["any"]]$true_minabs_causal_effect_mats[1:5],
    results[["any"]]$est_minabs_causal_effect_mats[1:5]
)
plotROCCurve(
    results[["maximal-arid"]]$true_minabs_causal_effect_mats[1:5],
    results[["maximal-arid"]]$est_minabs_causal_effect_mats[1:5]
)

options(vsc.dev.args = list(width = 250, height = 300))
pdf("./marg-tests/results/exp100any.pdf")
plotROCCurve(
    results100[["any"]]$true_minabs_causal_effect_mats,
    results100[["any"]]$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp100arid.pdf")
plotROCCurve(
    results100[["arid"]]$true_minabs_causal_effect_mats,
    results100[["arid"]]$est_minabs_causal_effect_mats
)
dev.off()
pdf("./marg-tests/results/exp100maximal-arid.pdf")
plotROCCurve(
    results100[["maximal-arid"]]$true_minabs_causal_effect_mats,
    results100[["maximal-arid"]]$est_minabs_causal_effect_mats
)
dev.off()



start <- proc.time()
results100projection <- run_experiment(N=100, aridities="projection", fast=TRUE)
end <-proc.time()
elapsed <- end - start
elapsed
saveRDS(results100projection, "./marg-tests/results/results100projection.RDS")

pdf("./marg-tests/results/exp100projection.pdf")
plotROCCurve(
    results100projection[["projection"]]$true_minabs_causal_effect_mats,
    results100projection[["projection"]]$est_minabs_causal_effect_mats
)
dev.off()

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
