% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjusted_envelopes.r
\name{GET.composite}
\alias{GET.composite}
\alias{dg.global_envelope_test}
\title{Adjusted global envelope tests}
\usage{
GET.composite(
  X,
  X.ls = NULL,
  nsim = 499,
  nsimsub = nsim,
  simfun = NULL,
  fitfun = NULL,
  calcfun = function(X) {
     X
 },
  testfuns = NULL,
  ...,
  type = "erl",
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater"),
  probs = c(0.025, 0.975),
  r_min = NULL,
  r_max = NULL,
  take_residual = FALSE,
  save.cons.envelope = savefuns,
  savefuns = FALSE,
  verbose = TRUE,
  MrkvickaEtal2017 = FALSE,
  mc.cores = 1L
)
}
\arguments{
\item{X}{An object containing the data in some form.
A \code{\link{curve_set}} object, or an \code{envelope} object of the \pkg{spatstat} package,
as the \code{curve_sets} argument of \code{\link{global_envelope_test}}
(need to provide \code{X.ls}), or
a fitted point process model of \pkg{spatstat} (e.g. object of class \code{ppm} or
\code{kppm}), or a point pattern object of class \code{ppp} of \pkg{spatstat},
or another data object (need to provide \code{simfun}, \code{fitfun}, \code{calcfun}).}

\item{X.ls}{A list of objects as \code{curve_sets} argument of \code{\link{global_envelope_test}},
giving the second stage simulations, see details.}

\item{nsim}{The number of simulations to be generated in the primary test.
Ignored if \code{X.ls} provided.}

\item{nsimsub}{Number of simulations in each basic test. There will be \code{nsim} repetitions
of the basic test, each involving \code{nsimsub} simulated realisations.
Total number of simulations will be nsim * (nsimsub + 1).}

\item{simfun}{A function for generating simulations from the null model. If given, this function
is called by \code{replicate(n=nsim, simfun(simfun.arg), simplify=FALSE)} to make nsim
simulations. Here \code{simfun.arg} is obtained by \code{fitfun(X)}.}

\item{fitfun}{A function for estimating the parameters of the null model.
The function should return the fitted model in the form that it can be directly
passed to \code{simfun} as its argument.}

\item{calcfun}{A function for calculating a summary function from a simulation of the model.
The default is the identity function, i.e. the simulations from the model are functions themselves.
The use of \code{calcfun} is still experimental. Preferably provide \code{X} and
\code{X.ls} instead, if \code{X} is not a point pattern or fitted point process model object
of \pkg{spatstat}.}

\item{testfuns}{A list of lists of parameters to be passed to the \code{envelope} function of \pkg{spatstat}
if \code{X} is a point pattern of a fitted point process model of \pkg{spatstat}.
A list of parameters should be provided for each test function that is to be used in the
combined test.}

\item{...}{Additional parameters to the \code{envelope} function of \pkg{spatstat} in the case where
only one test function is used. In that case, this is an alternative to providing the parameters in the
argument testfuns. If \code{envelope} is also used to generate simulations under the null
hypothesis (if simfun not provided), then also recall to specify how to generate the simulations.}

\item{type}{The type of the global envelope with current options for 'rank', 'erl', 'cont', 'area',
'qdir', 'st' and 'unscaled'. See details.}

\item{alpha}{The significance level. The 100(1-alpha)\% global envelope will be calculated
under the 'fwer' or 'fdr' control.
If a vector of values is provided, the global envelopes are calculated for each value.}

\item{alternative}{A character string specifying the alternative hypothesis.
Must be one of the following: "two.sided" (default), "less" or "greater".
The last two options only available for types \code{'rank'}, \code{'erl'},
\code{'cont'} and \code{'area'}.}

\item{probs}{A two-element vector containing the lower and upper
quantiles for the measure 'q' or 'qdir', in that order and on the interval [0, 1].
The default values are 0.025 and 0.975, suggested by Myllymäki et al. (2015, 2017).}

\item{r_min}{The minimum argument value to include in the test.}

\item{r_max}{The maximum argument value to include in the test.
in calculating functions by the \code{envelope} function of \pkg{spatstat}.}

\item{take_residual}{Logical. If TRUE (needed for visual reasons only) the mean of the simulated
functions is reduced from the functions in each first and second stage test.}

\item{save.cons.envelope}{Logical flag indicating whether to save the unadjusted envelope test results.}

\item{savefuns}{Logical flag indicating whether to save all the simulated function values.
Similar to the same argument of the \code{envelope} function of \pkg{spatstat}.}

\item{verbose}{Logical flag indicating whether to print progress reports during the simulations.
Similar to the same argument of \code{envelope} function of \pkg{spatstat}.}

\item{MrkvickaEtal2017}{Logical. If TRUE, type is "st" or "qdir" and several test functions are used,
then the combined scaled MAD envelope presented in Mrkvička et al. (2017) is calculated. Otherwise,
the two-step procedure described in \code{\link{global_envelope_test}} is used for combining the tests.
Default to FALSE. The option is kept for historical reasons.}

\item{mc.cores}{The number of cores to use, i.e. at most how many child processes will be run simultaneously.
Must be at least one, and parallelization requires at least two cores. On a Windows computer mc.cores must be 1
(no parallelization). For details, see \code{\link{mclapply}}, for which the argument is passed.
Parallelization can be used in generating simulations and in calculating the second stage tests.}
}
\value{
An object of class \code{global_envelope} or \code{combined_global_envelope}, which can be
printed and plotted directly. See \code{\link{global_envelope_test}}.
}
\description{
Adjusted global envelope tests for composite null hypothesis.
}
\details{
The specification of X, X.ls, fitfun, simfun is important:
\itemize{
\item If \code{X.ls} is provided, then the global envelope test is calculated based on
functions in these objects. \code{X} should be a \code{\link{curve_set}} object,
or an \code{envelope} object of \pkg{spatstat}, including the observed function and simulations
from the tested model. \code{X.ls} should be a list of \code{curve_set} or
envelope (of R package \pkg{spatstat}) objects, where each component contains an "observed"
function f that has been simulated from the model fitted to the data and the simulations
that have been obtained from the same model that has been fitted to the "observed" f.
The user has the responsibility that the functions have been generated correctly,
the test is done based on these provided simulations. See the examples.
\item Otherwise, if \code{simfun} and \code{fitfun} are provided, \code{X} can be general.
The function \code{fitfun} is used for fitting the desired model M and the function \code{simfun}
for simulating from a fitted model M. These functions should be coupled with each other such
that the object returned by \code{fitfun} is directly accepted as the (single) argument in
\code{simfun}.
In the case, that the global envelope should not be calculated directly for \code{X} (\code{X} is
not a function), \code{calcfun} can be used to specify how to calculate the function from
\code{X} and from simulations generated by \code{simfun}.
Special attention is needed in the correct specification of the functions, see examples.
 \item Otherwise, \code{X} should be either a fitted (point process) model object or a \code{ppp}
  object of the R package \pkg{spatstat}.
\itemize{
  \item If \code{X} is a fitted (point process) model object of the R package \pkg{spatstat},
then the simulations from this model are generated and summary functions for testing calculated
by the \code{envelope} function of \pkg{spatstat}. Which summary function to use and how to calculate it,
can be passed to \code{envelope} either in \code{...} or \code{testfuns}.
Unless otherwise specified the default function of \code{envelope},
i.g. the K-function, is used. The argument \code{testfuns} should be used to specify the
test functions in the case where one wants to base the test on several test functions.
  \item If \code{X} is a \code{ppp} object of \pkg{spatstat}, then the \code{envelope} function
is used for simulations and model fitting and the complete spatial randomness (CSR) is tested
(with intensity parameter).
}
}

For the rank envelope test, the global envelope test is the test described in
Myllymäki et al. (2017) with the adjustment of Baddeley et al. (2017).
For other test types, the test (also) uses the two-stage procedure of Dao and Genton (2014) with
the adjustment of Baddeley et al. (2017) as descripbed in Myllymäki and Mrkvička (2024).

See examples also in \code{\link{saplings}}.
}
\examples{
# Graphical normality test (Myllymaki and Mrkvicka, 2020, Section 3.3.)
#=========================
if(require("fda.usc", quietly=TRUE)) {
  data("poblenou")
  dat <- poblenou[['nox']][['data']][,'H10']
  n <- length(dat)

  # The number of simulations
  \donttest{nsim <- nsimsub <- 199}
  \dontshow{nsim <- nsimsub <- 19}

  set.seed(200127)
  # General setup
  #==============
  # 1. Fit the model
  mu <- mean(dat)
  sigma <- sd(dat)
  # 2. Simulate a sample from the fitted null model and
  #    compute the test vectors for data (obs) and each simulation (sim)
  r <- seq(min(dat), max(dat), length=100)
  obs <- stats::ecdf(dat)(r)
  sim <- sapply(1:nsimsub, function(i) {
    x <- rnorm(n, mean = mu, sd = sigma)
    stats::ecdf(x)(r)
  })
  cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))

  # 3. Simulate another sample from the fitted null model.
  # 4. Fit the null model to each of the patterns,
  #    simulate a sample from the null model,
  #    and compute the test vectors for all.
  cset.ls <- vector("list", nsim)
  for(rep in 1:nsim) {
    x <- rnorm(n, mean = mu, sd = sigma)
    mu2 <- mean(x)
    sigma2 <- sd(x)
    obs2 <- stats::ecdf(x)(r)
    sim2 <- sapply(1:nsimsub, function(i) {
      x2 <- rnorm(n, mean = mu2, sd = sigma2)
      stats::ecdf(x2)(r)
    })
    cset.ls[[rep]] <- create_curve_set(list(r = r, obs = obs2, sim_m = sim2))
  }
  # Perform the adjusted test
  res <- GET.composite(X = cset, X.ls = cset.ls, type = 'erl')
  plot(res) + ggplot2::labs(x = "NOx", y = "Ecdf")
}

# Example of a point pattern data
#================================
# Test the fit of a Matern cluster process.
\donttest{
if(require("spatstat.model", quietly=TRUE)) {
  data("saplings")
  saplings <- as.ppp(saplings, W = square(75))

  # First choose the r-distances
  rmin <- 0.3; rmax <- 10; rstep <- (rmax-rmin)/500
  r <- seq(0, rmax, by = rstep)

  # Number of simulations
  nsim <- 19 # Increase nsim for serious analysis!

  # Option 1: Give the fitted model object to GET.composite
  #--------------------------------------------------------
  # This can be done and is preferable when the model is
  # a point process model of spatstat.
  # 1. Fit the Matern cluster process to the pattern
  # (using minimum contrast estimation with the K-function)
  M1 <- kppm(saplings~1, clusters = "MatClust", statistic = "K")
  summary(M1)
  # 2. Make the adjusted global area rank envelope test using the L(r)-r function
  adjenvL <- GET.composite(X = M1, nsim = nsim,
              testfuns = list(L =list(fun="Lest", correction="translate",
                            transform=expression(.-r), r=r)), # passed to envelope
              type = "area", r_min = rmin, r_max = rmax)
  # Plot the test result
  plot(adjenvL)

  # Option 2: Generate the simulations "by yourself"
  #-------------------------------------------------
  # and provide them as curve_set or envelope objects
  # Preferable when you want to have a control
  # on the simulations yourself.
  # 1. Fit the model
  M1 <- kppm(saplings~1, clusters = "MatClust", statistic = "K")
  # 2. Generate nsim simulations by the given function using the fitted model
  X <- envelope(M1, nsim = nsim, savefuns = TRUE,
                fun = "Lest", correction = "translate",
                transform = expression(.-r), r = r)
  plot(X)
  # 3. Create another set of simulations to be used to estimate
  # the second-state p-value (as proposed by Baddeley et al., 2017).
  simpatterns2 <- simulate(M1, nsim = nsim)
  # 4. Calculate the functions for each pattern
  simf <- function(rep) {
    # Fit the model to the simulated pattern Xsims[[rep]]
    sim_fit <- kppm(simpatterns2[[rep]], clusters = "MatClust", statistic = "K")
    # Generate nsim simulations from the fitted model
    envelope(sim_fit, nsim = nsim, savefuns = TRUE,
             fun = "Lest", correction = "translate",
             transform = expression(.-r), r = r)
  }
  X.ls <- parallel::mclapply(X = 1:nsim, FUN = simf, mc.cores = 1) # list of envelope objects
  # 5. Perform the adjusted test
  res <- GET.composite(X = X, X.ls = X.ls, type = "area",
                      r_min = rmin, r_max = rmax)
  plot(res)
}}
}
\references{
Baddeley, A., Hardegen, A., Lawrence, T., Milne, R. K., Nair, G. and Rakshit, S. (2017). On two-stage Monte Carlo tests of composite hypotheses. Computational Statistics and Data Analysis 114: 75-87. doi: http://dx.doi.org/10.1016/j.csda.2017.04.003

Dao, N.A. and Genton, M. (2014). A Monte Carlo adjusted goodness-of-fit test for parametric models describing spatial point patterns. Journal of Graphical and Computational Statistics 23, 497-517.

Mrkvička, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27(5): 1239-1255. doi: 10.1007/s11222-016-9683-9

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381-404. doi: 10.1111/rssb.12172

Myllymäki, M. and Mrkvička, T. (2024). GET: Global envelopes in R. Journal of Statistical Software 111(3), 1-40. doi: 10.18637/jss.v111.i03
}
\seealso{
\code{\link{global_envelope_test}}, \code{\link{plot.global_envelope}}, \code{\link{saplings}}
}
