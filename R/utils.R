compute_SBC_cache_blocks <- function(datasets, ..., block_size, cache_prefix, cache_suffix = ".rds") {
  starts <- seq(1, length(datasets), by = block_size)
  ends <- c(starts[2:length(starts)] - 1, length(datasets))
  res_list <- list()
  for(i in 1:length(starts)) {
    res_list[[i]] <- compute_SBC(datasets[starts[i]:ends[i]], ..., cache_mode = "results", cache_location = paste0(cache_prefix, i, cache_suffix))
  }
  do.call(bind_results, res_list)
}

get_partial_cache_blocks <- function(datasets, backend, ..., block_size, cache_prefix, cache_suffix = ".rds") {
  starts <- seq(1, length(datasets), by = block_size)
  ends <- c(starts[2:length(starts)] - 1, length(datasets))
  res_list <- list()
  backend_hash <- SBC_backend_hash_for_cache(backend)
  for(i in 1:length(starts)) {
    filename <- paste0(cache_prefix, i, cache_suffix)
    if(file.exists(filename)) {
      results_from_cache <- readRDS(filename)
      data_hash <- rlang::hash(datasets[starts[i]:ends[i]])

      # TODO create SBC::is_cache_valid
      result <- NULL
      if(!is.list(results_from_cache) ||
         !all(
           c("result", "backend_hash", "data_hash", "thin_ranks", "dquants","keep_fits")
           %in% names(results_from_cache))) {
        warning("Cache file ", filename, " exists but is in invalid format. Ignoring.")
      } else {
        if(results_from_cache$backend_hash != backend_hash) {
          message("Cache file ", filename, " exists but the backend hash differs. Use with caution!")
        } else if(results_from_cache$data_hash != data_hash) {
        message("Cache file exists ", filename, " but the datasets hash differs. Use with caution!")
        }
        result <- tryCatch(validate_SBC_results(results_from_cache$result),
                           error = function(e) { NULL })
      }
      if(!is.null(result)) {
        res_list[[length(res_list) + 1]] <- result
      }
    }
  }
  do.call(bind_results, res_list)

}
