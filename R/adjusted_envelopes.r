# A helper function to generate simulated functions from "X" using given functions/arguments,
# when X is a ppp or model object of spatstat.
#
# a) If X is a point pattern (object of class "ppp"), then CSR is simulated (with variable number of points).
# b) If X is a fitted point process model (object of class "ppm" or "kppm"),
# the simulated model is that model (by spatstat).
# If testfuns is given, then also several different test functions can be calculated.
funcs_from_X_and_funs <- function(X, nsim, testfuns=NULL, ...,
                                  savepatterns=FALSE, verbose=FALSE, calc_funcs=TRUE) {
  extraargs <- list(...)
  nfuns <- length(testfuns)
  if(nfuns < 1) nfuns <- 1
  for(i in 1:nfuns)
    if(any(names(extraargs) %in% names(testfuns[[i]])))
      stop(paste("formal argument(s) \"", names(extraargs)[which(names(extraargs) %in% names(testfuns[[i]]))], "\" matched by multiple actual arguments.", sep=""))

  X.ls <- NULL
  # Simulations
  # Calculate the first test functions and generate simulations
  X.ls[[1]] <- do.call(spatstat::envelope, c(list(Y=X, nsim=nsim, simulate=NULL),
                          testfuns[[1]],
                          list(savefuns=TRUE, savepatterns=savepatterns, verbose=verbose),
                          extraargs))
  # More than one test function, calculate the rest if calc_funcs==TRUE
  if(calc_funcs & nfuns > 1) {
    simpatterns <- attr(X.ls[[1]], "simpatterns")
    for(i in 2:nfuns) {
      X.ls[[i]] <- do.call(spatstat::envelope, c(list(Y=X, nsim=nsim, simulate=simpatterns),
                              testfuns[[i]],
                              list(savefuns=TRUE, savepatterns=FALSE, verbose=verbose),
                              extraargs))
    }
  }
  X.ls
}

# A helper function to perform simulations for the GET.composite
#' @importFrom stats update
#' @importFrom parallel mclapply
adj.simulate <- function(X, nsim = 499, nsimsub = nsim,
                         simfun=NULL, fitfun=NULL, calcfun=function(X) { X },
                         testfuns=NULL, ..., verbose=TRUE, mc.cores=1L) {
  # Case 1: fitfun, simfun, calcfun provided
  # Model fitted by fitfun, simulated by simfun; X can be general
  if(!is.null(fitfun) & !is.null(simfun)) {
    message("Note: Model to be fitted by fitfun(X), simulations by simfun and calcfun;\n",
        "simfun should accept the object returned by fitfun as its argument.\n",
        "calcfun should accept the object returned by simfun as its argument.")
    # Fit the model to X
    simfun.arg <- fitfun(X) # fitted model to be passed to simfun
    # Generate nsim simulations by the given function using the fitted model
    stage1_cset_ls <- replicate(n=nsimsub, expr=simfun(simfun.arg), simplify=FALSE) # list of data objects
    stage1_cset_ls <- sapply(stage1_cset_ls, FUN=calcfun, simplify=TRUE) # matrix of functions
    stage1_cset_ls <- list(create_curve_set(list(r=1:nrow(stage1_cset_ls), obs=calcfun(X), sim_m=stage1_cset_ls)))
    # Create another set of simulations to be used to estimate the second-state p-value
    # following Baddeley et al. (2017).
    Xsims <- replicate(n=nsim, expr=simfun(simfun.arg), simplify=FALSE) # list of data objects
    loopfun <- function(rep) {
      # Fit the model to the simulated pattern Xsims[[rep]]
      simfun.arg <- fitfun(Xsims[[rep]])
      # Generate nsimsub simulations from the fitted model
      cset <- replicate(n=nsimsub, expr=simfun(simfun.arg), simplify=FALSE) # list of data objects
      cset <- sapply(cset, FUN=calcfun, simplify=TRUE) # matrix of functions
      list(create_curve_set(list(r=1:nrow(cset), obs=calcfun(X), sim_m=cset)))
    }
    stage2_csets_lsls <- parallel::mclapply(X=1:nsim, FUN=loopfun, mc.cores=mc.cores)
  }
  # Case 2: simfun or fitfun not provided; X should be a ppp or model object of spatstat
  # a) If X is a ppp object, the tested model is CSR
  # b) If X is a model object of spatstat, then spatstat is used for fitting and simulation.
  else {
    # Check if X is a (ppp) model object of spatstat
    Xispppmodel <- spatstat::is.ppm(X) || spatstat::is.kppm(X) || spatstat::is.lppm(X) || spatstat::is.slrm(X)
    if(!spatstat::is.ppp(X) & !Xispppmodel) stop("fitfun or simfun not provided and X is not a ppp nor a fitted model object of spatstat.")
    if(Xispppmodel) message("X is a fitted model object of spatstat;\n using spatstat to simulate the model and calculate the test functions.")
    else
      message("Note: \'simfun\' and/or \'fitfun\' not provided and \'X\' is a ppp object of spatstat.\n",
          "The spatstat's function \'envelope\' is used for simulations and model fitting, \n",
          "and CSR is tested (with intensity parameter).")
    # Create simulated functions from the given model
    stage1_cset_ls <- funcs_from_X_and_funs(X, nsim=nsimsub, testfuns=testfuns, ...,
                                            savepatterns=FALSE, verbose=verbose, calc_funcs=TRUE)
    # Create another set of simulations to be used to estimate the second-state p-value
    # following Baddeley et al. (2017).
    simpatterns <- attr(funcs_from_X_and_funs(X, nsim=nsim, testfuns=testfuns, ...,
                                              savepatterns=TRUE, verbose=verbose, calc_funcs=FALSE)[[1]],
                        "simpatterns")
    # For each of the simulated patterns in 'simpatterns', fit the model and
    # create nsim simulations from it
    loopfun <- function(rep) {
      # Create a simulation
      Xsim <- simpatterns[[rep]]
      if(Xispppmodel) Xsim <- update(X, Xsim)
      # Create simulations from the given model and calculate the test functions
      funcs_from_X_and_funs(Xsim, nsim=nsimsub, testfuns=testfuns, ...,
                            savepatterns=FALSE, verbose=verbose, calc_funcs=TRUE)
    }
    stage2_csets_lsls <- mclapply(X=1:nsim, FUN=loopfun, mc.cores=mc.cores)
  }
  list(X=stage1_cset_ls, X.ls=stage2_csets_lsls)
}

GE.attr <- function(x, name = "p", ...) {
  if(inherits(x, c("global_envelope", "global_envelope_2d"))) a <- attr(x, name)
  if(inherits(x, c("combined_global_envelope", "combined_global_envelope_2d"))) a <- attr(attr(x, "level2_ge"), name)
  a
}

# A helper function to make a global envelope test for the purposes of GET.composite
# @param curve_sets_ls A list of curve_sets.
# @inheritParams GET.composite
adj.GET_helper <- function(curve_sets, type, alpha, alternative, ties, probs, MrkvickaEtal2017, ..., save.envelope=FALSE) {
  if(length(curve_sets) > 1 & MrkvickaEtal2017 & type %in% c("st", "qdir")) {
    global_envtest <- combined_scaled_MAD_envelope_test(curve_sets, type=type, alpha=alpha, probs=probs, ...)
  }
  else {
    global_envtest <- global_envelope_test(curve_sets, type=type, alpha=alpha,
                                           alternative=alternative, ties=ties, probs=probs, ...)
  }
  res <- structure(list(stat = as.numeric(GE.attr(global_envtest, "M")[1]), pval = GE.attr(global_envtest, "p")),
                   class = "adj_GET_helper")
  if(save.envelope) attr(res, "envelope_test") <- global_envtest
  res
}

#' Adjusted global envelope tests
#'
#' Adjusted global envelope tests for composite null hypothesis.
#'
#'
#' The specification of X, X.ls, fitfun, simfun is important:
#' \itemize{
#' \item If \code{X.ls} is provided, then the global envelope test is calculated based on
#' functions in these objects. \code{X} should be a \code{curve_set} (see \code{\link{create_curve_set}})
#' or an \code{envelope} object of \pkg{spatstat} including the observed function and simulations
#' from the tested model. \code{X.ls} should be a list of \code{curve_set} or
#' envelope (of R library \pkg{spatstat}) objects, where each component contains an "observed"
#' function f that has been simulated from the model fitted to the data and the simulations
#' that have been obtained from the same model that has been fitted to the "observed" f.
#' The user has the responsibility that the functions have been generated correctly,
#' the test is done based on these provided simulations. See the examples.
#' \item Otherwise, if \code{simfun} and \code{fitfun} are provided, \code{X} can be general.
#' The function \code{fitfun} is used for fitting the desired model M and the function \code{simfun}
#' for simulating from a fitted model M. These functions should be coupled with each other such
#' that the object returned by \code{fitfun} is directly accepted as the (single) argument in
#' \code{simfun}.
#' In the case, that the global envelope should not be calculated directly for \code{X} (\code{X} is
#' not a function), \code{calcfun} can be used to specify how to calculate the function from
#' \code{X} and from simulations generated by \code{simfun}.
#' Special attention is needed in the correct specification of the functions, see examples.
#'  \item Otherwise, \code{X} should be either a fitted (point process) model object or a \code{ppp}
#'   object of the R library \pkg{spatstat}.
#' \itemize{
#'   \item If \code{X} is a fitted (point process) model object of the R library \pkg{spatstat},
#' then the simulations from this model are generated and summary functions for testing calculated
#' by the \code{envelope} function of \pkg{spatstat}. Which summary function to use and how to calculate it,
#' can be passed to \code{envelope} either in \code{...} or \code{testfuns}.
#' Unless otherwise specified the default function of \code{envelope},
#' i.g. the K-function, is used. The argument \code{testfuns} should be used to specify the
#' test functions in the case where one wants to base the test on several test functions.
#'   \item If \code{X} is a \code{ppp} object of \pkg{spatstat}, then the \code{envelope} function
#' is used for simulations and model fitting and the complete spatial randomness (CSR) is tested
#' (with intensity parameter).
#' }
#' }
#'
#' For the rank envelope test, the global envelope test is the test described in
#' Myllymäki et al. (2017) with the adjustment of Baddeley et al. (2017).
#' For other test types, the test (also) uses the two-stage procedure of Dao and Genton (2014) with
#' the adjustment of Baddeley et al. (2017) as descripbed in Myllymäki and Mrkvička (2020).
#'
#' See examples also in \code{\link{saplings}}.
#'
#' @param X An object containing the data in some form.
#' A \code{curve_set} (see \code{\link{create_curve_set}}) or an \code{envelope}
#' object (of the \pkg{spatstat} package), as the \code{curve_sets} argument of \code{\link{global_envelope_test}}
#' (need to provide \code{X.ls}), or
#' a fitted point process model of \pkg{spatstat} (e.g. object of class \code{ppm} or
#' \code{kppm}), or a point pattern object of class \code{ppp} of \pkg{spatstat},
#' or another data object (need to provide \code{simfun}, \code{fitfun}, \code{calcfun}).
#' @param X.ls A list of objects as \code{curve_sets} argument of \code{\link{global_envelope_test}},
#' giving the second stage simulations, see details.
#' @param nsim The number of simulations to be generated in the primary test.
#' Ignored if \code{X.ls} provided.
#' @param nsimsub Number of simulations in each basic test. There will be \code{nsim} repetitions
#' of the basic test, each involving \code{nsimsub} simulated realisations.
#' Total number of simulations will be nsim * (nsimsub + 1).
#' @param simfun A function for generating simulations from the null model. If given, this function
#' is called by \code{replicate(n=nsim, simfun(simfun.arg), simplify=FALSE)} to make nsim
#' simulations. Here \code{simfun.arg} is obtained by \code{fitfun(X)}.
#' @param fitfun A function for estimating the parameters of the null model.
#' The function should return the fitted model in the form that it can be directly
#' passed to \code{simfun} as its argument.
#' @param calcfun A function for calculating a summary function from a simulation of the model.
#' The default is the identity function, i.e. the simulations from the model are functions themselves.
#' The use of \code{calcfun} is still experimental. Preferably provide \code{X} and
#' \code{X.ls} instead, if \code{X} is not a point pattern or fitted point process model object
#' of \pkg{spatstat}.
#' @param testfuns A list of lists of parameters to be passed to the \code{envelope} function of \pkg{spatstat}
#' if \code{X} is a point pattern of a fitted point process model of \pkg{spatstat}.
#' A list of parameters should be provided for each test function that is to be used in the
#' combined test.
#' @param ... Additional parameters to the \code{envelope} function of \pkg{spatstat} in the case where
#' only one test function is used. In that case, this is an alternative to providing the parameters in the
#' argument testfuns. If \code{envelope} is also used to generate simulations under the null
#' hypothesis (if simfun not provided), then also recall to specify how to generate the simulations.
#' @inheritParams global_envelope_test
#' @param r_min The minimum argument value to include in the test.
#' @param r_max The maximum argument value to include in the test.
#' in calculating functions by the \code{envelope} function of \pkg{spatstat}.
#' @param take_residual Logical. If TRUE (needed for visual reasons only) the mean of the simulated
#' functions is reduced from the functions in each first and second stage test.
#' @param save.cons.envelope Logical flag indicating whether to save the unadjusted envelope test results.
#' @param savefuns Logical flag indicating whether to save all the simulated function values.
#' Similar to the same argument of the \code{envelope} function of \pkg{spatstat}.
#' @param verbose Logical flag indicating whether to print progress reports during the simulations.
#' Similar to the same argument of \code{envelope} function of \pkg{spatstat}.
#' @param MrkvickaEtal2017 Logical. If TRUE, type is "st" or "qdir" and several test functions are used,
#' then the combined scaled MAD envelope presented in Mrkvička et al. (2017) is calculated. Otherwise,
#' the two-step procedure described in \code{\link{global_envelope_test}} is used for combining the tests.
#' Default to FALSE. The option is kept for historical reasons.
#' @param mc.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously.
#' Must be at least one, and parallelization requires at least two cores. On a Windows computer mc.cores must be 1
#' (no parallelization). For details, see \code{\link{mclapply}}, for which the argument is passed.
#' Parallelization can be used in generating simulations and in calculating the second stage tests.
#' @return An object of class \code{global_envelope} or \code{combined_global_envelope}, which can be
#' printed and plotted directly. See \code{\link{global_envelope_test}}.
#' @references
#' Baddeley, A., Hardegen, A., Lawrence, T., Milne, R. K., Nair, G. and Rakshit, S. (2017). On two-stage Monte Carlo tests of composite hypotheses. Computational Statistics and Data Analysis 114: 75-87. doi: http://dx.doi.org/10.1016/j.csda.2017.04.003
#'
#' Dao, N.A. and Genton, M. (2014). A Monte Carlo adjusted goodness-of-fit test for parametric models describing spatial point patterns. Journal of Graphical and Computational Statistics 23, 497-517.
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27(5): 1239-1255. DOI: 10.1007/s11222-016-9683-9
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381-404. doi: 10.1111/rssb.12172
#'
#' Myllymäki, M. and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#'
#' @seealso \code{\link{global_envelope_test}}, \code{\link{plot.global_envelope}}, \code{\link{saplings}}
#' @export
#' @aliases dg.global_envelope_test
#' @importFrom parallel mclapply
#' @importFrom stats quantile
#' @examples
#' # Graphical normality test (Myllymaki and Mrkvicka, 2020, Section 3.3.)
#' #=========================
#' if(require("fda.usc", quietly=TRUE)) {
#'   data("poblenou")
#'   dat <- poblenou[['nox']][['data']][,'H10']
#'   n <- length(dat)
#'
#'   # The number of simulations
#'   \donttest{nsim <- nsimsub <- 199}
#'   \dontshow{nsim <- nsimsub <- 19}
#'
#'   set.seed(200127)
#'   # General setup
#'   #==============
#'   # 1. Fit the model
#'   mu <- mean(dat)
#'   sigma <- sd(dat)
#'   # 2. Simulate a sample from the fitted null model and
#'   #    compute the test vectors for data (obs) and each simulation (sim)
#'   r <- seq(min(dat), max(dat), length=100)
#'   obs <- stats::ecdf(dat)(r)
#'   sim <- sapply(1:nsimsub, function(i) {
#'     x <- rnorm(n, mean = mu, sd = sigma)
#'     stats::ecdf(x)(r)
#'   })
#'   cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))
#'
#'   # 3. Simulate another sample from the fitted null model.
#'   # 4. Fit the null model to each of the patterns,
#'   #    simulate a sample from the null model,
#'   #    and compute the test vectors for all.
#'   cset.ls <- list()
#'   for(rep in 1:nsim) {
#'     x <- rnorm(n, mean = mu, sd = sigma)
#'     mu2 <- mean(x)
#'     sigma2 <- sd(x)
#'     obs2 <- stats::ecdf(x)(r)
#'     sim2 <- sapply(1:nsimsub, function(i) {
#'       x2 <- rnorm(n, mean = mu2, sd = sigma2)
#'       stats::ecdf(x2)(r)
#'     })
#'     cset.ls[[rep]] <- create_curve_set(list(r = r, obs = obs2, sim_m = sim2))
#'   }
#'   # Perform the adjusted test
#'   res <- GET.composite(X=cset, X.ls=cset.ls, type='erl')
#'   plot(res) + ggplot2::labs(x="NOx", y="Ecdf")
#' }
#'
#' @examples
#' # Example of a point pattern data
#' #================================
#' # Test the fit of a Matern cluster process.
#' \donttest{
#' if(require("spatstat", quietly=TRUE)) {
#'   data("saplings")
#'   saplings <- as.ppp(saplings, W=square(75))
#'
#'   # First choose the r-distances
#'   rmin <- 0.3; rmax <- 10; rstep <- (rmax-rmin)/500
#'   r <- seq(0, rmax, by=rstep)
#'
#'   # Number of simulations
#'   nsim <- 19 # Increase nsim for serious analysis!
#'
#'   # Option 1: Give the fitted model object to GET.composite
#'   #--------------------------------------------------------
#'   # This can be done and is preferable when the model is
#'   # a point process model of spatstat.
#'   # 1. Fit the Matern cluster process to the pattern
#'   # (using minimum contrast estimation with the K-function)
#'   M1 <- kppm(saplings~1, clusters = "MatClust", statistic="K")
#'   summary(M1)
#'   # 2. Make the adjusted global area rank envelope test using the L(r)-r function
#'   adjenvL <- GET.composite(X = M1, nsim = nsim,
#'               testfuns = list(L = list(fun="Lest", correction="translate",
#'                              transform = expression(.-r), r=r)), # passed to envelope
#'               type = "area", r_min=rmin, r_max=rmax)
#'   # Plot the test result
#'   plot(adjenvL)
#'
#'   # Option 2: Generate the simulations "by yourself"
#'   #-------------------------------------------------
#'   # and provide them as curve_set or envelope objects
#'   # Preferable when you want to have a control
#'   # on the simulations yourself.
#'   # 1. Fit the model
#'   M1 <- kppm(saplings~1, clusters = "MatClust", statistic="K")
#'   # 2. Generate nsim simulations by the given function using the fitted model
#'   X <- envelope(M1, nsim=nsim, savefuns=TRUE,
#'                 fun="Lest", correction="translate",
#'                 transform = expression(.-r), r=r)
#'   plot(X)
#'   # 3. Create another set of simulations to be used to estimate
#'   # the second-state p-value (as proposed by Baddeley et al., 2017).
#'   simpatterns2 <- simulate(M1, nsim=nsim)
#'   # 4. Calculate the functions for each pattern
#'   simf <- function(rep) {
#'     # Fit the model to the simulated pattern Xsims[[rep]]
#'     sim_fit <- kppm(simpatterns2[[rep]], clusters = "MatClust", statistic="K")
#'     # Generate nsim simulations from the fitted model
#'     envelope(sim_fit, nsim=nsim, savefuns=TRUE,
#'              fun="Lest", correction="translate",
#'              transform = expression(.-r), r=r)
#'   }
#'   X.ls <- parallel::mclapply(X=1:nsim, FUN=simf, mc.cores=1) # list of envelope objects
#'   # 5. Perform the adjusted test
#'   res <- GET.composite(X=X, X.ls=X.ls, type="area",
#'                       r_min=rmin, r_max=rmax)
#'   plot(res)
#' }}
GET.composite <- function(X, X.ls = NULL,
                         nsim = 499, nsimsub = nsim,
                         simfun=NULL, fitfun=NULL, calcfun=function(X) { X },
                         testfuns=NULL, ...,
                         type = "erl",
                         alpha = 0.05, alternative = c("two.sided","less", "greater"),
                         probs = c(0.025, 0.975),
                         r_min=NULL, r_max=NULL, take_residual=FALSE,
                         save.cons.envelope = savefuns, savefuns = FALSE,
                         verbose=TRUE, MrkvickaEtal2017 = FALSE, mc.cores=1L) {
  alt <- match.arg(alternative)

  # Simulations
  #------------
  # 1) All simulations provided
  #----------------------------
  if(!is.null(X.ls)) {
    # Check dimensions of each curve_set and transform to curve_sets
    if(!(length(class(X)) == 1 && class(X) == "list")) {
      X <- list(X) # The observed curves (a list of curve_set objects)
      X.ls <- lapply(X.ls, FUN=list) # The simulated curves (a list of lists of curve_set objects)
    }
    picked_attr_ls <- lapply(X, FUN=pick_attributes, alternative=alt)
    X <- check_curve_set_dimensions(X)
    X.ls <- lapply(X.ls, FUN = check_curve_set_dimensions)
    # Check equality of dimensions over repetitions
    if(!all(sapply(X.ls, FUN=function(curve_set) { identical(curve_set[[1]]$r, y=X[[1]]$r) }))) stop("The number of argument values in the observed and simulated sets of curves differ.")
    if(!all(sapply(X.ls, FUN=function(curve_set) { curve_set_nfunc(curve_set[[1]]) == curve_set_nfunc(X[[1]]) }))) stop("The number of simulations in the observed and simulated sets of curves differ.")
    # Checking r_min, r_max
    if(!is.null(r_min) & length(r_min) != length(X)) stop("r_min should give the minimum distances for each of the test functions.")
    if(!is.null(r_max) & length(r_max) != length(X)) stop("r_max should give the maximum distances for each of the test functions.")
    message("Using the simulations provided in X and X.ls.")
  }
  # 2) Simulations if X.ls not provided
  #------------------------------------
  else { # Perform simulations
    if(verbose) message("Performing simulations, ...")
    tmp <- adj.simulate(X=X, nsim=nsim, nsimsub=nsimsub,
                       simfun=simfun, fitfun=fitfun, calcfun=calcfun, testfuns=testfuns, ...,
                       verbose=verbose, mc.cores=mc.cores)
    picked_attr_ls <- lapply(tmp$X, FUN=pick_attributes, alternative=alt)
    X <- check_curve_set_dimensions(tmp$X)
    X.ls <- lapply(tmp$X.ls, FUN = check_curve_set_dimensions)
  }

  # Crop curves (if r_min or r_max given)
  #------------
  nfuns <- length(X)
  # Cropping (and tranforming to curve_sets)
  if(!is.null(r_min) | !is.null(r_max)) {
    tmpfunc <- function(i, csets) crop_curves(csets[[i]], r_min=r_min[i], r_max=r_max[i])
    X <- lapply(1:nfuns, FUN=tmpfunc, csets=X)
    for(sim in 1:length(X.ls)) X.ls[[sim]] <- lapply(1:nfuns, FUN=tmpfunc, csets=X.ls[[sim]])
  }

  # Take residual
  #--------------
  if(take_residual) {
    tmpfunc <- function(i, csets) residual(csets[[i]])
    X <- lapply(1:nfuns, FUN=tmpfunc, csets=X)
    for(sim in 1:length(X.ls)) X.ls[[sim]] <- lapply(1:nfuns, FUN=tmpfunc, csets=X.ls[[sim]])
  }

  # Individual tests
  #-----------------
  # For data
  obs_res <- adj.GET_helper(curve_sets=X, type=type, alpha=alpha, alternative=alt, ties="midrank",
                           probs=probs, MrkvickaEtal2017=MrkvickaEtal2017, save.envelope=TRUE)
  # For simulations
  loopfun <- function(rep) {
    tmp <- adj.GET_helper(curve_sets=X.ls[[rep]], type=type, alpha=alpha, alternative=alt, ties="midrank",
                         probs=probs, MrkvickaEtal2017=MrkvickaEtal2017)
    list(stat = tmp$stat, pval = tmp$pval)
  }
  mclapply_res <- mclapply(X=1:length(X.ls), FUN=loopfun, mc.cores=mc.cores)
  stats <- sapply(mclapply_res, function(x) x$stat)
  pvals <- sapply(mclapply_res, function(x) x$pval)

  # The adjusted test
  #------------------
  if(MrkvickaEtal2017 & type %in% c("st", "qdir") & nfuns > 1) { # Combined tests as in Mrkvicka et al. (2017)
    res <- attr(obs_res, "envelope_test")
    #-- The rank test at the second level
    # Calculate the critical rank / alpha
    Malpha_star <- quantile(stats, probs=alpha, type=1)
    all_curves <- data_and_sim_curves(attr(res, "level2_curve_set")) # the second step "curves"
    nr <- ncol(all_curves)
    Nfunc <- nrow(all_curves)
    LB <- array(0, nr)
    UB <- array(0, nr)
    for(i in 1:nr){
      Hod <- sort(all_curves[,i])
      LB[i]<- Hod[Malpha_star]
      UB[i]<- Hod[Nfunc-Malpha_star+1]
    }
    # -> Malpha_star, LB, UB of the (second level) rank test
    # Update res object with adjusted values
    attr(res, "level2_ge")$lo <- LB
    attr(res, "level2_ge")$hi <- UB
    attr(attr(res, "level2_ge"), "M_alpha_star") <- Malpha_star # Add Malpha_star
    # Re-calculate the new qdir/st envelopes
    envchars <- combined_scaled_MAD_bounding_curves_chars(X, type = type)
    central_curves_ls <- lapply(X, function(x) get_T_0(x))
    bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=UB,
                                                           lower_f=envchars$lower_f, upper_f=envchars$upper_f)
    # Update the first level envelopes for plotting
    for(i in 1:length(res)) {
      res[[i]]$lo <- bounding_curves$lower_ls[[i]]
      res[[i]]$hi <- bounding_curves$upper_ls[[i]]
    }
  }
  else { # Otherwise, the ERL test is used at the second level of a combined test
    if(type == "rank" & nfuns == 1) {
      # Calculate the critical rank (instead of alpha) and the adjusted envelope following Myllymäki et al. (2017)
      Malpha_star <- quantile(stats, probs=alpha, type=1)
      all_curves <- data_and_sim_curves(X[[1]]) # all the functions
      nr <- length(X[[1]]$r)
      Nfunc <- nrow(all_curves)
      LB <- array(0, nr)
      UB <- array(0, nr)
      for(i in 1:nr){
        Hod <- sort(all_curves[,i])
        LB[i]<- Hod[Malpha_star]
        UB[i]<- Hod[Nfunc-Malpha_star+1]
      }
      # For getting automatically an global_envelope object, first call central_region
      res <- central_region(X, type=type, coverage=1-alpha, alternative=alt, central="mean")
      # Update res object with adjusted values
      res$lo <- LB
      res$hi <- UB
      attr(res, "M_alpha_star") <- Malpha_star # Add Malpha_star
    }
    else {
      alpha_star <- quantile(pvals, probs=alpha, type=1)
      res <- global_envelope_test(X, type=type, alpha=alpha_star, alternative=alt)
      # Original p-value
      attr(res, "p_original") <- attr(res, "p")
      # Adjusted p-value
      attr(res, "p") <- estimate_p_value(x=-obs_res$pval, sim_vec=-pvals, ties="conservative")
      # Save additional arguments
      if(nfuns == 1) {
        attr(res, "alpha") <- alpha
        attr(res, "alpha_star") <- alpha_star
      }
      else {
        attr(res, "alpha") <- attr(attr(res, "level2_ge"), "alpha") <- alpha
        attr(res, "alpha_star") <- attr(attr(res, "level2_ge"), "alpha_star") <- alpha_star
      }
    }
  }

  # Update attributes
  if(nfuns == 1) {
    for(n in names(picked_attr_ls[[1]])) attr(res, n) <- picked_attr_ls[[1]][[n]]
  }
  else { # nfuns > 1
    for(n in names(picked_attr_ls[[1]]))
      for(i in 1:nfuns)
        attr(res[[i]], n) <- picked_attr_ls[[i]][[n]]
  }
  attr(res, "method") <- "Adjusted global test" # Change method name
  attr(res, "call") <- match.call() # Update "call" attribute
  # Additions
  if(savefuns) {
    attr(res, "simfuns") <- X
    attr(res, "simfuns.ls") <- X.ls
  }
  if(save.cons.envelope) attr(res, "unadjusted_test") <- attr(obs_res, "envelope_test")
  attr(res, "simulated_p") <- pvals
  attr(res, "simulated_k") <- stats
  # Return
  res
}
