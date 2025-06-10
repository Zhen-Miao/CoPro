
## ---------- soft-deprecated wrappers for backward compatibility ----------
#' @noRd
newCoPro <- function(...) {
  .Deprecated("newCoProSingle")           # nudge users
  newCoProSingle(...)
}

#' @noRd
newCoProm <- function(...) {
  .Deprecated("newCoProMulti")
  newCoProMulti(...)
}


#' subsetDataOne
#' @param ... Arguments to be passed to the new function
#' @returns The CoPro object
#' @noRd
subsetDataOne <- function(...) {
  .Deprecated("subsetData")
  subsetData(...)
}


#' @param ... Arguments to be passed to the new function
#' @returns The CoPro object
#' @noRd
subsetDataMulti <- function(...) {
  .Deprecated("subsetData")
  subsetData(...)
}

#' @noRd
computePCAMulti <- function(...) {
  .Deprecated("computePCA")
  computePCA(...)
}
