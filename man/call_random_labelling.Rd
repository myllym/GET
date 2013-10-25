\name{call_random_labelling}
\alias{call_random_labelling}
\title{Helper function for both \code{\link{random_labelling_K_f}} and
\code{\link{random_labelling_L_f}}.}
\usage{
  call_random_labelling(pattern, do_besags_L,
    mtf_name = "m", r_max = NULL, r_vec = NULL,
    nsim = 999L, calc_theo = TRUE,
    edge_correction = "translate", method = "permute",
    use_biased_lambda2 = FALSE)
}
\arguments{
  \item{calc_theo}{Whether to calculate the theoretical
  summary function as well.}

  \item{do_besags_L}{Whether to calculate the K_f or the
  L_f function.}

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

  \item{nsim}{The number of permutations.}

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
}
\value{
  A list with either two or three components. 'obs' has the
  summary function that was asked for. 'r' contains the
  radius values. 'unmarked' contains the unmarked K or L
  function, if asked for.
}
\description{
  Helper function for both
  \code{\link{random_labelling_K_f}} and
  \code{\link{random_labelling_L_f}}.
}

