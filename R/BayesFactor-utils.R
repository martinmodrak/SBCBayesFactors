#' A helper function that collapses variables corresponding to arrays in the `lmBF` output for neater visuals.
combine_lmBF_arrays <- function(x) {
  indices_removed <- gsub("-[^-]*$", "", x)
  unique_arrays <- sort(unique(indices_removed))
  res <- list()
  for(i in 1:length(unique_arrays)) {
    res[[unique_arrays[i]]] <- x[indices_removed == unique_arrays[i]]
  }
  res
}
