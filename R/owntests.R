source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/api.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/greedy.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/icfmagCustom.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/plotting.R", encoding = "UTF-8")
source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/simulation.R", encoding = "UTF-8")

source("/home/sebastiaan/Documents/Masterscriptie/greedyBAPs/R/margs.R", encoding = "UTF-8")

set.seed(1)
bap <- generateUniformBAP(p = 4, N = 1)[[1]]
params <- generateBAPParameters(bap)
data <- generateData(n = 100, params = params)
plotBAP(bap, nodesize=40)
pairs(data)

ancestor_bap <- t(cbind(
  c(0, 1, 100, 0, 100),
  c(0, 0, 1, 100, 0),
  c(0, 0, 0, 1, 0),
  c(0, 0, 0, 0, 0),
  c(0, 1, 0, 0, 0)
))
plotBAP(ancestor_bap)
ancestors(ancestor_bap, 4) # should be 4 3 2 1 5
children(ancestor_bap, 5)
parents(ancestor_bap, 2) # should be 1 5
descendants(ancestor_bap, 5) # should be 5 2 3 4
siblings(ancestor_bap, 1) # should be 3 5
siblings(ancestor_bap, 4) # should be 2
district(ancestor_bap, 3) # should be 3 1 5

reachable_closure(ancestor_bap, c(3,5))
reachable_closure(ancestor_bap, 4)

closure_bap <- t(cbind(
    c(0, 1, 100),
    c(0, 0, 1),
    c(0, 0, 0)
))
plotBAP(closure_bap)
reachable_closure(closure_bap, c(2, 3)) # should be 1 2 3, cause 1 isn't fixable
reachable_closure(closure_bap, 1)
reachable_closure(closure_bap, 2)
reachable_closure(closure_bap, 3)
reachable_closure(closure_bap, c(1,3))

maximal_arid_projection(closure_bap)

# weird as hell, used for hashing, I think
connectedComponents(ancestor_bap)$comp

marg_ex1 <- t(cbind(
    c(0, 100, 1, 100),
    c(0, 0, 100, 1),
    c(0, 0, 0, 0),
    c(0, 0, 0, 0)
))
plotBAP(marg_ex1)
plotBAP(maximal_arid_projection(marg_ex1))
maximal_arid_projection(marg_ex1)


# testing a graph that I think is bow-free, but not arid
bowfree_also_arid <- t(cbind(
    c(0, 1, 100, 0),
    c(0, 0, 1, 100),
    c(0, 0, 0, 1),
    c(0, 0, 0, 0)
))
plotBAP(bowfree_also_arid)
lapply(c(1,2,3,4), reachable_closure, graph=bowfree_also_arid)

bowfree_not_arid <- t(cbind(
    c(0, 1, 100, 0),
    c(0, 0, 0, 100),
    c(0, 1, 0, 100),
    c(1, 0, 0, 0)
))
bowfree_not_arid <- t(cbind(
    c(0, 1, 100, 100),
    c(0, 0, 1, 100),
    c(0, 0, 0, 0),
    c(0, 0, 1, 0)
))
plotBAP(bowfree_not_arid)
lapply(c(1,2,3,4), reachable_closure, graph=bowfree_not_arid)

interesting_marg <- maximal_arid_projection(bowfree_not_arid)
plotBAP(interesting_marg)

ancestors_not_arid <- lapply(1:4, ancestors, graph=bowfree_not_arid)
closures_not_arid <- lapply(1:4, reachable_closure, graph=bowfree_not_arid)
reachable_closure(bowfree_not_arid, c(2,4))
district(bowfree_not_arid, c(2))
parents(bowfree_not_arid, 1) # closure of 1
parents(bowfree_not_arid, 1:4)
