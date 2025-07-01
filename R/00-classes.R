#' @import Matrix
#' @import methods
# Define necessary class unions from the original CoPro object
setClassUnion("matrixOrSparseMatrix", c("matrix", "dgCMatrix", "dgTMatrix"))
setClassUnion("factorOrCharacter", c("factor", "character"))
setClassUnion("matrixOrDataFrame", c("matrix", "data.frame"))

# Define a virtual class that is a union of 'matrix' and 'sparseMatrix'
#' CoPro object of spatial transcriptomics data
#' @import Matrix
#'
#' @slot normalizedData A `matrix` object to store normalized data.
#' @slot normalizedDataSub A `matrix` object to store the subset of
#' the normalized data. The subset only contain relevant cell types
#' of interest specified by the user.
#' @slot locationData A `data.frame` object to store the location. It should
#' either contain two columns named by "x" and "y", or three columns named by
#' "x", "y", and "z". No other names allowed
#' @slot locationDataSub A `data.frame` object to store the subset of
#' the location data.
#' @slot metaData A `data.frame` object to store metadata for each cell.
#' @slot metaDataSub A `data.frame` object to store the subset of
#' the meta data.
#' @slot cellTypes A `vector` object with elements being character. It should
#' match the number of cells in the data matrix and each represents a cell type
#' label of a cell.
#' @slot cellTypesSub A `vector` object with elements being character. It stores
#' the subset of  the cell type labels.
#' @slot cellTypesOfInterest A `vector` object with elements being character.
#' Specifies the cell types of interest. It will be used to subset the dataset
#' @slot pcaResults A `list` object storing PCA results after integration.
#'  Recommended structure: `list(slideID = list(cellType = pc_matrix))`.
#' @slot pcaGlobal A `list` object storing PCA results for each cell type.
#' @slot distances A `list` object to store the pairwise distances between any
#' two cell types of interest.
#' @slot geneList A `vector` object with elements being character. To store the
#' gene names.
#' @slot kernelMatrices A `list` object. To store the kernel matrix generated
#' from the distance matrices.
#' @slot sigmaValues A `vector` object with elements being numeric. To store
#' a set of sigma values used for generating the kernel matrix.
#' @slot nPCA A single numeric value. Number of PCs to retain for downstream
#' analyses.
#' @slot nCC A single numeric value. Number of canonical components to retain
#' for downstream analyses.
#' @slot scalePCs A `logical` value. Whether to scale each PC before computing
#'  skrCCA
#' @slot skrCCAOut A `list` object. Output from the skrCCA.
#' @slot skrCCAPermuOut A `list` object. Output from the skrCCA after
#' permutation. This helps establish the null distribution
#' @slot cellPermu A `list` object that stores the cell permutation labels
#' @slot nPermu A `numeric` value specifying the number of permutations
#'  conducted.
#' @slot cellScores A `matrix` object. Cell scores for each cell type.
#' @slot geneScores A `matrix` object. Gene scores for each cell type.
#' @slot geneScoresTest A `list` object. Tested gene scores
#' @slot normalizedCorrelation A `list` object. Normalized correlation values
#' for each sigma value.
#' @slot normalizedCorrelationPermu A `list` object.
#'  Normalized correlation values
#' for each sigma value after permutation
#' @slot sigmaValueChoice A `numeric` value. The optimal sigma squared based
#' on the median normalized correlation value.
#' @name CoPro-class
setClass("CoPro",
         contains  = "VIRTUAL",
         slots = list(

           ## cell by gene data matrix
           normalizedData = "matrixOrSparseMatrix",
           normalizedDataSub = "matrixOrSparseMatrix",
           integratedData = "list", # Stores output of integrateSlidesMulti

           ## location data
           locationData = "matrixOrDataFrame",
           locationDataSub = "matrixOrDataFrame",

           ## cell by meta.data matrix
           metaData = "data.frame",
           metaDataSub = "data.frame",

           ## cell type labels for each cell
           cellTypes = "factorOrCharacter",
           cellTypesSub = "factorOrCharacter",

           ## cell types of interest
           cellTypesOfInterest = "character",

           ## pca results of normalized data matrix
           pcaResults = "list",
           pcaGlobal = "list",      # Stores list(cellType = pc_result)
           ## pcaGlobal is only used in multi-slide setting
           ## pcaResults are also defined differently  --> let us fix this
           ## In single slide setting, let us just use pcaGlobal for consistency
           ## there is no need to have pcaResults in single slide setting

           ## geneList
           geneList = "character",

           ## some parameters
           nPCA = "numeric",
           scalePCs = "logical",

           # Distance & Kernel (Slide-Specific if multiple slideID)
           distances = "list",      # Stores list(slideID = list(ct1 = list(ct2 = dist)))
           kernelMatrices = "list", # Stores list(sigma = list(slideID = list(ct1 = list(ct2 = K))))
           sigmaValues = "numeric",

           # skrCCA Output (Shared Weights, Slide-Specific Scores)
           nCC = "numeric",
           skrCCAOut = "list",      # Stores list(sigma = list(cellType = W)) - Shared
           normalizedCorrelation = "list", # Store per-slide or aggregate results
           cellScores = "list",     # Stores list(sigma = list(slideID = list(cellType = Scores)))
           geneScores = "list",     # Stores list(sigma = list(cellType = gene_scores)) - Potentially Shared
           geneScoresTest = "list",
           sigmaValueChoice = "numeric",

           ## permutation output
           nPermu = "numeric",
           skrCCAPermuOut = "list",
           cellPermu = "list",
           normalizedCorrelationPermu = "list"
         ),
         prototype = list(
           cellTypesOfInterest = character(0),
           geneList = character(0),
           sigmaValues = numeric(0),
           nPCA = numeric(0),
           nCC = numeric(0),
           nPermu = numeric(0),
           sigmaValueChoice = numeric(0),
           scalePCs = logical(0)
         )
)


setClass(
  "CoProSingle",
  contains = "CoPro"
)


setClass(
  "CoProMulti",
  contains = "CoPro",
  slots = list(
    slideList = "character"
  ),
  prototype = list(
    slideList = character(0)
  )
)

setClass("CoProm", contains = "CoProMulti")    # alias for *old* CoProm

