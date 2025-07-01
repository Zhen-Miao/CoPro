
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
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  # Check if RcppAnnoy is available
  if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
    stop(paste("Package 'RcppAnnoy' is needed for this function to work.",
               "Please install it."), call. = FALSE)
  }

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

# Search an Annoy approximate nearest neighbor index
#
# @param Annoy index, built with AnnoyBuildIndex
# @param query A set of data to be queried against the index
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time trade off between better accuracy and speed.
# @param include.distance Include the corresponding distances in the result
#
# @return A list with 'nn.idx' (for each element in 'query', the index of the
# nearest k elements in the index) and 'nn.dists' (the distances of the nearest
# k elements)
#
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#
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


########################## modified ##########################
#' Create an Annoy index
#'
#' @note Function exists because it's not exported from \pkg{uwot}
#'
#' @param name Distance metric name
#' @param ndim Number of dimensions
#'
#' @return An nn index object
#'
CreateAnn <- function(name, ndim) {
  # Check if RcppAnnoy is available
  if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
    stop(paste("Package 'RcppAnnoy' is needed for this function to work.",
         "Please install it."), call. = FALSE)
  }

  return(switch(
    EXPR = name,
    cosine = methods::new(Class = RcppAnnoy::AnnoyAngular, ndim),
    manhattan = methods::new(Class = RcppAnnoy::AnnoyManhattan, ndim),
    euclidean = methods::new(Class = RcppAnnoy::AnnoyEuclidean, ndim),
    hamming = methods::new(Class = RcppAnnoy::AnnoyHamming, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  ))
}


########################## modified ##########################
#' Internal helper function to dispatch to various neighbor finding methods
#'
#' @param data Input data
#' @param query Data to query against data
#' @param k Number of nearest neighbors to compute
#' @param method Nearest neighbor method to use: "rann", "annoy"
#' @param cache.index Store algorithm index with results for reuse
#' @param ... additional parameters to specific neighbor finding method
#
NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  # Check if SeuratObject is available
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop( paste("Package 'SeuratObject' is needed for this function to work.",
         "Please install it."), call. = FALSE)
  }

  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
  results <- do.call(what = 'AnnoyNN', args = args)

  # Use getClass to access the class definition from SeuratObject
  n.ob <- methods::new(
    Class = methods::getClass("Neighbor", where = asNamespace("SeuratObject")),
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )

  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    methods::slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
}


#' Main function to compute KNN matrix
#' @param data Matrix of data to query against object.
#' @param query query
#' @param search.k search.k
#' @param k number of neighbors
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean,
#' cosine, manhattan, and hamming
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param nn.eps Error bound when performing nearest neighbor search using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param cache.index Include cached index in returned Neighbor object
#' (only relevant if return.neighbor = TRUE)
#' @param index Precomputed index. Useful if querying new data against existing
#' index to avoid recomputing.
#'
#' @importFrom RANN nn2
#' @importFrom methods as
#'
#' @export
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

  # Check if SeuratObject is available
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop( paste("Package 'SeuratObject' is needed for this function to work.",
                "Please install it."), call. = FALSE)
  }

  nn.ranked <- SeuratObject::Indices(object = nn.ranked)
  # convert nn.ranked into a Graph
  j <- as.numeric(x = t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), nrow(x = object))), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = object)
  colnames(x = nn.matrix) <- rownames(x = object)

  return(nn.matrix)
}

