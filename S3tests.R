pheromone_at.pheromoneMatrix <- function(mat, row, col) {
    return(mat[row, col])
}
pheromone_at.pheromoneList <- function(lst, row, col) {
    if (paste(row, col) %in% names(lst)) {
        return(lst[[paste(row, col)]])
    } else {
        return(0)
    }
}
pheromone_at <- function(container, row, col) {
    UseMethod(generic="pheromone_at", object=container)#, row, col)
}


pheromone_mat <- diag(1, nrow=3, ncol=4)
pheromone_list <- list("1 1"=2, "2 2"=2, "3 3"=2, "1 4"=5)
# this doesn't work:
# pheromone_list$class <- "pheromoneList"

class(pheromone_mat) <- "pheromoneMatrix"
class(pheromone_list) <- "pheromoneList"

pheromone_to_mat <- function(container, nrow, ncol) {
    result <- matrix(0, nrow, ncol)
    for (irow in 1:nrow) {
        for (icol in 1:ncol) {
            # print(paste(irow, icol))
            result[irow, icol] <- pheromone_at(container, irow, icol)
        }
    }
    return(result)
}

pheromone_to_mat(pheromone_mat, 3, 4)
pheromone_to_mat(pheromone_list, 3, 4)

