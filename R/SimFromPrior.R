# SimFromPrior
# Daniel J. Schad, Dez. 13, 2018

SimFromPrior <- function(priors,class="b",coef="") {
  priors_par <- priors$prior[priors$class==class & priors$coef==coef]
  priors_par <- strsplit(priors_par,"(",fixed=TRUE)[[1]]
  par <- strsplit(priors_par[2],")")[[1]]
  par <- as.numeric(strsplit(par,",")[[1]])
  priors_par <- priors_par[1]
  if (priors_par=="normal") {
    #cat(paste0("rnorm(1,",par[1],",",par[2],")\n"))
    parSamp <- rnorm(1,par[1],par[2])
    if (class%in%c("sigma","sd"))
      while (parSamp<=0) parSamp <- rnorm(1,par[1],par[2])
  }
  if (priors_par=="lkj") {
    library(rethinking)
    #cat(paste0("rethinking::rlkjcorr(1,2,",par[1],")\n"))
    parSamp <- rethinking::rlkjcorr(1,2,par[1])[1,2]
  }
  parSamp
}

