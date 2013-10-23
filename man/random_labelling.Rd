\name{random_labelling}
\alias{random_labelling}
\title{Random labeling with K_f functions}
\usage{
  random_labelling(pattern, mtf_name = "m", nsim = 999L,
    r_min = NULL, r_max = NULL, r_vec = NULL, use_L = TRUE,
    edge_correction = "translate", use_theo = TRUE,
    method = "permute", use_biased_lambda2 = FALSE)
}
\arguments{
  \item{pattern}{A \code{\link[spatstat]{ppp}} object as
  the simple marked point pattern to be analysed. The marks
  need to be in the form of a numeric vector. The window
  has to have the type "rectangle".}

  \item{mtf_name}{A vector of mark test function names. "1"
  stands for the unmarked K-function. Accepted values are
  '1', 'm', 'mm', 'gamma', 'gammaAbs' and 'morAbs'.}

  \item{r_max}{A positive scalar value representing the
  maximum radius that should be considered. r_vec overrides
  r_max. By default, r_max is NULL and will get a sensible
  default.}

  \item{r_vec}{A monotonically increasing vector of
  non-negative r-values to act as the endpoints of the bins
  for the K_f-functions. r_vec overrides r_max. The bins
  are exclusive on the left and inclusive on the right. If
  the first vector element has value zero, it will be
  regarded as the collapsed bin [0, 0], and the next bin
  will start from and exclude 0.}

  \item{edge_correction}{The name of the edge correction to
  be used. Options are 'translate' and 'none'.}

  \item{method}{The name of the method to create
  simulations under the null hypothesis. 'permute' results
  in permutations of the marks. Using 'sample' will sample
  the marks from the empirical mark distribution with
  replacement. 'permute' is the default.}

  \item{use_biased_lambda2}{A logical scalar on whether to
  use the biased or the unbiased (in the Poisson case)
  estimate of the intensity squared.}

  \item{r_min}{The minimum radius to include.}

  \item{use_theo}{Whether to use the theoretical summary
  function or the mean of the simulations.}
}
\description{
  Perform simulations under the random labelling for a
  marked point pattern and estimate K_f- or L_f-functions
  for the pattern and simulations.
}
\references{
  Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D.
  (2013). Deviation test construction and power comparison
  for marked spatial point patterns. arXiv:1306.1028
}
\seealso{
  \code{\link{rank_envelope}}, \code{\link{st_envelope}},
  \code{\link{qdir_envelope}}, \code{\link{deviation_test}}
}

