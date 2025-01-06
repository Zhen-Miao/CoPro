
## Acknowledgment:  the following functions were modified based off Seurat
## function for finding KNN.

# Run annoy
#
# @param data Data to build the index with
# @param query A set of data to be queried against data
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @param include.distance Include the corresponding distances
# @param index optional index object, will be recomputed if not provided
#
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL
) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}

# Build the annoy index
#
# @param data Data to build the index with
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
#
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
#
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}


AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = plan(), what = "multicore")) {
    oplan <- plan(strategy = "sequential")
    on.exit(plan(oplan), add = TRUE)
  }
  res <- future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}


# Create an Annoy index
#
# @note Function exists because it's not exported from \pkg{uwot}
#
# @param name Distance metric name
# @param ndim Number of dimensions
#
# @return An nn index object
#
#' @importFrom methods new
#' @importFrom RcppAnnoy AnnoyAngular AnnoyManhattan AnnoyEuclidean AnnoyHamming
#
CreateAnn <- function(name, ndim) {
  return(switch(
    EXPR = name,
    cosine = new(Class = AnnoyAngular, ndim),
    manhattan = new(Class = AnnoyManhattan, ndim),
    euclidean = new(Class = AnnoyEuclidean, ndim),
    hamming = new(Class = AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
}

# Internal helper function to dispatch to various neighbor finding methods
#
# @param data Input data
# @param query Data to query against data
# @param k Number of nearest neighbors to compute
# @param method Nearest neighbor method to use: "rann", "annoy"
# @param cache.index Store algorithm index with results for reuse
# @param ... additional parameters to specific neighbor finding method
#
#' @importFrom methods new
#' @importClassesFrom SeuratObject Neighbor
#
NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
  results <- do.call(what = 'AnnoyNN', args = args)

  n.ob <- new(
    Class = 'Neighbor',
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
}



# Main function to compute KNN matrix
find_neighbor_mat <- function(data, k = 20, nn.method = "annoy", query = NULL,
                           annoy.metric = "euclidean", nn.eps = 0,
                           cache.index = FALSE, index = NULL,
                           n.trees = 50, search.k = NULL) {

  object <- as.matrix(data)
  query <- query %||% object
  k.param <- k

  ## require row names
  if(is.null(rownames(data))){
    stop("Row names must be given and match cell names.")
  }

  ## compute nn.ranked
  nn.ranked <- NNHelper(
    data = object,
    query = query,
    k = k,
    method = nn.method,
    n.trees = n.trees,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric,
    cache.index = cache.index,
    index = index
  )

  nn.ranked <- Indices(object = nn.ranked)
  # convert nn.ranked into a Graph
  j <- as.numeric(x = t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), nrow(x = object))), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = object)
  colnames(x = nn.matrix) <- rownames(x = object)

  return(nn.matrix)
}

