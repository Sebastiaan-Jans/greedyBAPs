source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/api.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/greedy.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/icfmagCustom.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/plotting.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/simulation.R", encoding = "UTF-8")

source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/margs.R", encoding = "UTF-8")
library(foreach) # for %dopar% (paralellism)
library(ggm)

set.seed(1)
bap <- generateUniformBAP(p = 4, N = 1)[[1]]
params <- generateBAPParameters(bap)
data <- generateData(n = 10000, params = params)
# plotBAP(bap)
# pairs(data)
score_graph(bap, data)

res <- greedySearch(cov(data), nrow(data), n.restarts=10, mc.cores=1, aridity="arid")
plotScoreCurves(res$all.scores, res$times, true_graph_score=score_graph(bap, data))

marg_res <- greedySearch(cov(data), nrow(data), n.restarts=10,
    # verbose=TRUE,
    mc.cores=4,
    margs.only=TRUE)
options(vsc.dev.args = list(width = 600, height = 600))
plotScoreCurves(marg_res$all.scores, marg_res$times, true_graph_score=score_graph(bap, data))

# normal non-marg ROC results
effects_res <- causalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10
)
roc <- plotROCCurve(effects_res$CE.gt, effects_res$CE.greedy)

effects_res2 <- causalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10
)
roc <- plotROCCurve(effects_res2$CE.gt, effects_res2$CE.greedy)

effects_marg_res <- causalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10,
    margs.only = TRUE
)
roc <- plotROCCurve(effects_marg_res$CE.gt, effects_marg_res$CE.greedy)


## Custom simulations:

# normal non-marg ROC results
effects_res <- customCausalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10
)
roc <- plotROCCurve(effects_res$CE.gt, effects_res$CE.greedy)

effects_res2 <- customCausalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10
)
roc <- plotROCCurve(effects_res2$CE.gt, effects_res2$CE.greedy)

effects_marg_res <- customCausalEffectsSimulation(
    N = 10,
    p = 5, # number of variables, I think
    n = 1000,
    mc.cores = 4,
    n.restarts = 10,
    max.in.degree = 2,
    max.steps = 100,
    max.iter.ricf = 10,
    equivalent.eps = 1e-10,
    margs.only = TRUE
)
roc <- plotROCCurve(effects_marg_res$CE.gt, effects_marg_res$CE.greedy)

# give a somewhat functional stack trace
options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  if (!interactive()) {
    q(status = 1)
  }
})

