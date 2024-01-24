node_names <- function(graph) {
    if (!is.null(colnames(graph))) return(colnames(graph))
    return(1:ncol(graph))
}

# Return a version of the graph where the bi-edges are stored in both transpose index positions
# example: the edge 2 <-> 4 implies graph[2, 4] == graph[4, 2] == 100
make_symmetric <- function(graph) {
    biedge_mask <- graph == 100
    symmetric_mask <- pmax((biedge_mask), t(biedge_mask))
    return(100 * symmetric_mask + (graph * (graph != 100)))
}

transpose_100 <- function(matrix) {
  matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(matrix)]
  return(matrix)
}


# just work with integer vertices for now
maximal_arid_projection <- function(graph) {
    nodes <- 1:ncol(graph)

    # initialize a new graph with the same vertices but no edges
    new_graph <- graph
    new_graph[] <- 0

    # TODO: make this an actual dict of reachable closures
    # format: v=reachable_closure(v)
    reachable_closures = list()
    ancestors <- lapply(nodes, ancestors, graph=graph)
    for (a in nodes) {
        for (b in tail(nodes, -a)) {
            u <- NULL
            v <- NULL
            closure <- NULL

            if (a %in% ancestors[[b]]) {
                u <- a
                v <- b
            } else if (b %in% ancestors[[a]]) {
                u <- b
                v <- a
            }
            
            added_edge <- FALSE
            if (!is.null(u)) {
                if (!(v %in% names(reachable_closures))) {
                    reachable_closures[[v]] <- reachable_closure(graph, v)$closure
                }
                closure <- reachable_closures[[v]]
                if (u %in% parents(graph, closure)) {
                    new_graph[u, v] = 1
                    added_edge <- TRUE
                }
            }
            if (!added_edge) { # semantically superfluous check, but saves time
                combined_closure_result <- reachable_closure(graph, c(a, b))
                combined_closure <- combined_closure_result$closure

                b_district <- district(graph, b)
                # if (a %in% b_district && is_subset(combined_closure, b_district)) { % wrong!
                if (is_bidirected_connected(graph, combined_closure)) {
                    new_graph[a, b] <- new_graph[b, a] <- 100
                }
            }
        }
    }
    return(new_graph)
}

# helper function that returns whether vector smaller is contained within larger
# note that it also returns TRUE when smaller equals larger
is_subset <- function(smaller, larger) {
    all(smaller %in% larger)
}

reachable_closure <- function(graph, vertices, already_fixed=c(), get_names = FALSE) {
    # we'll deal with integer id variables for now
    fixed_vertices <- already_fixed
    graph_vertices <- 1:ncol(graph)

    remaining_vertices <- setdiff(
        graph_vertices,
        union(vertices, already_fixed)
    )
    fixing_order <- c()
    fixed <- TRUE
    new_graph <- graph

    while (length(remaining_vertices > 0) && fixed) {
        fixed <- FALSE
        for (vert in remaining_vertices) {
            if (length(
                intersect(descendants(new_graph, vert), district(new_graph, vert))
            ) == 1) {
                fix_result <- fix(new_graph, vert, already_fixed=fixed_vertices)
                new_graph <- fix_result$graph
                fixed_vertices <- fix_result$fixed_vertices
                remaining_vertices <- setdiff(remaining_vertices, vert)
                fixing_order <- c(fixing_order, vert)
                fixed <- TRUE
                break
            }
        }
    }

    reachable_closure <- setdiff(graph_vertices, fixed_vertices)
    if (get_names) {
        reachable_closure <- node_names(graph)[reachable_closure]
        fixing_order <- node_names(graph)[fixing_order]
    }
    return(list(closure=reachable_closure, fixing_order=fixing_order, graph=new_graph))
}

fix <- function(graph, vertex, already_fixed=c()) {
    fixed_vertices <- c(already_fixed, vertex)
    # remove all incoming directed edges (and possible bidirected)
    graph[, vertex] <- 0
    # remove any remaining bidirected edges from the row
    graph[vertex, ][siblings(graph, vertex)] <- 0
    # note that ananke actually stores the districts in the object, and recomputes here
    return(list(graph=graph, fixed_vertices=fixed_vertices))
}


# works for any number of variables by taking the union over each's parents
parents <- function(graph, variables, get_names = FALSE) {
    parents_onevar <- function(var) which(graph[, var] == 1)
    pars <- Reduce(union, sapply(variables, parents_onevar))
    if (!get_names) return(pars)
    return(node_names(graph)[pars])
}

# works for any number of variables by taking the union over each's parents
children <- function(graph, variables, get_names = FALSE) {
    children_onevar <- function(var) which(graph[var, ] == 1)
    children <- Reduce(union, sapply(variables, children_onevar))
    if (!get_names) return(children)
    return(node_names(graph)[children])
}

siblings <- function(graph, variables, get_names = FALSE) {
    # bi-edges are (often?) only marked once, so we need to check the column *and* the row
    siblings_onevar <- function(var) {
        return(unique(c(
            which(graph[, var] == 100),
            which(graph[var, ] == 100)
        )))
    }
    siblings <- Reduce(union, lapply(variables, siblings_onevar))
    if (!get_names) return(siblings)
    return(node_names(graph)[siblings])
}

ancestors <- function(graph, variables, get_names = FALSE) {
    pars <- parents(graph, variables)
    ancs <- unique(c(variables, pars, unlist(lapply(pars, ancestors, graph=graph))))
    if (!get_names) return(ancs)
    return(node_names(graph)[ancs])
}

descendants <- function(graph, variables, get_names = FALSE) {
    childs <- children(graph, variables)
    descs <- unique(c(variables, childs, unlist(lapply(childs, descendants, graph=graph))))
    if (!get_names) return(descs)
    return(node_names(graph)[descs])
}

# bing actually improved my version a bit, this one visits the same state less
district <- function(graph, variable, get_names = FALSE) {
    district <- variable
    handled <- c(variable)

    district_inner <- function(variable) {
        for (neighbor in siblings(graph, variable)) {
            if (!(neighbor %in% handled)) {
                district <<- c(district, neighbor)
                handled <<- c(handled, neighbor)
                district_inner(neighbor)
            }
        }
    }

    district_inner(variable)
    district <- unique(district)

    if (!get_names) return(district)
    return(node_names(graph)[district])
}

# Return the subgraph containing only the vars in `variables`,
# leaving colnames and rownames intact
subgraph <- function(graph, variables) {
    sub <- graph[variables, variables]
    colnames(sub) <- rownames(sub) <- variables
    return(sub)
}

# Return whether a set of variables is bidirected-connected in graph
# if safe is TRUE, do not assume bi-edges are correctly stored twice, but costs a copy
# as.character is required because the subgraph matrix maintains the
# original indices, which is impossible with non-string indices
is_bidirected_connected <- function(graph, variable_set, safe=TRUE) {
    sub <- subgraph(graph, variable_set)
    if (safe) {
        sub <- make_symmetric(sub)
    }
    # The subgraph must only contain a single district,
    # checked by taking the district of one and seeing if it contains all
    # index <- as.character(which(colnames(graph) == variable_set[1])

    return(length(district(sub, 1)) == length(variable_set))
}
