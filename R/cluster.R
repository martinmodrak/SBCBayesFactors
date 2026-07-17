setup_cluster <- function(max_local_workers = Inf, max_server_workers = 140, force = FALSE) {
  is_local <- parallel::detectCores() < 80 # detectCores counts even unavailable cores
  if(is_local) {
    n_workers <- min(parallelly::availableCores(), max_local_workers)
    # future::plan(multisession, workers = n_workers)
  } else {
    n_workers <- min(parallelly::availableCores(), max_server_workers)
    # mirai::daemons(n_workers, dispatcher = FALSE, force = force)
    # future::plan(future.mirai::mirai_cluster)
  }
  withr::with_envvar(new = c("RENV_CONFIG_STARTUP_QUIET" = "true", "RENV_CONFIG_SANDBOX_ENABLED" = "FALSE"),
    {
      #if(is_local) {
        future::plan(future::multisession, workers = n_workers)
      #} else {
      #  mirai::daemons(n_workers, dispatcher = TRUE, force = force)
      #  future::plan(future.mirai::mirai_cluster)
      #}
    }
                     )
  options(SBC.generator_chunk_size = 500)
}
