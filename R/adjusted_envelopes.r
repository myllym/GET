# A helper function to generate simulations from "X" and calculate test functions for them.
#
# If X is a point pattern (object of class "ppp"), then CSR is simulated.
# If X is a fitted point process model (object of class "ppm" or "kppm"),
# the simulated model is that model.
# simfun and simfun.arg can alternatively be used to specify how to generate the simulations,
# see global_envelope_with_sims below.
# If testfuns is given, then also several different test functions can be calculated.
simpatterns_and_funcs_from_X <- function(X, nsim, simfun=NULL, simfun.arg=NULL, testfuns=NULL, ...,
                                         savepatterns=FALSE, verbose=FALSE, calc_funcs=TRUE) {
  nfuns <- length(testfuns)
  testfuns_argnames_ls <- lapply(testfuns, function(x) names(x))
  X.ls <- NULL
  if(!is.null(simfun)) {
    # Create simulations by the given function
    simpatterns <- replicate(n=nsim, simfun(simfun.arg), simplify=FALSE)
    # Calculate the test functions
    if(is.null(testfuns)) {
      X.ls[[1]] <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ..., savefuns = TRUE, savepatterns = savepatterns, verbose=verbose)
    }
    else {
      # More than one test function, calculate only the first one unless if calc_funcs==TRUE
      i <- 1
      tmpstring <- paste("X.ls[[1]] <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ", sep="")
      for(j in 1:length(testfuns_argnames_ls[[i]]))
        tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], "=testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
      tmpstring <- paste(tmpstring, "savefuns=TRUE, savepatterns=savepatterns, verbose=verbose, ...)", sep="")
      eval(parse(text = tmpstring))
      if(calc_funcs & nfuns > 1) {
        for(i in 2:nfuns) {
          tmpstring <- paste("X.ls[[i]] <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ", sep="")
          for(j in 1:length(testfuns_argnames_ls[[i]]))
            tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], "=testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
          tmpstring <- paste(tmpstring, "savefuns=TRUE, savepatterns=FALSE, verbose=verbose, ...)", sep="")
          eval(parse(text = tmpstring))
        }
      }
    }
  }
  else {
    # Create simulations from the given model and calculate the test functions
    if(is.null(testfuns)) {
      X.ls[[1]] <- spatstat::envelope(X, nsim=nsim, ..., savefuns = TRUE, savepatterns = savepatterns, verbose=verbose)
    }
    else {
      i <- 1
      tmpstring <- paste("X.ls[[1]] <- spatstat::envelope(X, nsim=nsim, ", sep="")
      for(j in 1:length(testfuns_argnames_ls[[i]]))
        tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], " = testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
      tmpstring <- paste(tmpstring, "savefuns = TRUE, savepatterns = TRUE, verbose=verbose, ...)", sep="")
      eval(parse(text = tmpstring))
      if(calc_funcs) {
        simpatterns <- attr(X.ls[[1]], "simpatterns")
        if(calc_funcs & nfuns > 1) {
          for(i in 2:nfuns) {
            tmpstring <- paste("X.ls[[i]] <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ", sep="")
            for(j in 1:length(testfuns_argnames_ls[[i]]))
              tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], "=testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
            tmpstring <- paste(tmpstring, "savefuns=TRUE, savepatterns=FALSE, verbose=verbose, ...)", sep="")
            eval(parse(text = tmpstring))
          }
        }
      }
    }
  }
  X.ls
}

# A global envelope test
#
# A global envelope test including simulations from a point process model.
#
#
# Details of the tests are given in \code{\link{global_envelope_test}}.
# For the combined directional quantile and studentized envelope tests of Mrkvička et al. (2017)
# see \code{\link{combined_scaled_MAD_envelope}}.
#
# The specification of X is important here:
# if simfun is not provided, the function \code{\link[spatstat]{envelope}} is used to generate
# simulations under the null hypothesis and to calculate the test functions (specified in the
# arguments ...) and then
# \itemize{
# \item If X is a point pattern, the null hypothesis is CSR.
# \item If X is a fitted model, the null hypothesis is that model.
# }
# If simfun is provided, then the null model is the one simulated by this given function,
# and X is expected to be a point pattern of \code{\link[spatstat]{ppp}} object, for which data
# pattern and simulations \code{\link[spatstat]{envelope}} calculates the test statistics.
#
# If savefuns is TRUE, all the simulated functions are saved in an attribute "simfuns" as a list
# of curve sets for each test function.
#
# @param X An object containing point pattern data. A point pattern (object of class "ppp")
# or a fitted point process model (object of class "ppm" or "kppm"). See
# \code{\link[spatstat]{envelope}}.
# @param nsim The number of simulations.
# @param simfun A function for generating simulations from the null model. If given, this function
# is called by replicate(n=nsim, simfun(simfun.param), simplify=FALSE) to make nsim simulations.
# The function should return an \code{\link[spatstat]{ppp}} object as those are further passed to
# \code{\link[spatstat]{envelope}}.
# If the function is not provided, then \code{\link[spatstat]{envelope}} is used also for generating
# the point patterns from the null hypothesis.
# @param simfun.arg The parameter to be passed to simfun. The function simfun should handle
# with the structure of simfun.param.
# @param testfuns A list of lists of parameters to be passed to \code{\link[spatstat]{envelope}}.
# A list of parameters should be provided for each test function that is to be used in the combined test.
# @param ... Additional parameters to \code{\link[spatstat]{envelope}} in the case where only one test
# function is used. In that case, this is an alternative to providing the parameter in the argument
# testfuns.
# For example, the test function in the argument 'fun' and further specifications regarding that.
# If \code{\link[spatstat]{envelope}} is also used to generate simulations under the null hypothesis
# (if simfun not provided), then also recall to specify how to generate the simulations.
# @param type The type of the test to be applied, see \code{\link{global_envelope_test}}.
# @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
# @param alternative A character string specifying the alternative hypothesis. Must be one of
# the following: "two.sided" (default), "less" or "greater" for "rank". Relevant only for the
# rank test (otherwise ignored).
# @param r_min The minimum radius to include in the test.
# @param r_max The maximum radius to include in the test. Note: cannot be larger than r-values used
# in calculating functions by \code{\link[spatstat]{envelope}}.
# @param take_residual If (needed for visual reasons only) the theoretical or mean behaviour of the
# test function is reduced from the test functions. If TRUE, then: If \code{\link[spatstat]{envelope}}
# provides the theoretical value 'theo' for the simulated model, then this value is used. Otherwise,
# the theoretical function is taken as the mean of the simulations.
# @param ties Ties method to be passed to \code{\link{global_envelope_test}}.
# Used to obtain a point estimate for the p-value for the "rank" test. Default to extreme rank
# length p-value.
# @param save.envelope Logical flag indicating whether to save all envelope test results.
# @param savefuns Logical flag indicating whether to save all the simulated function values.
# See \code{\link[spatstat]{envelope}}.
# @param savepatterns Logical flag indicating whether to save all the simulated point patterns.
# See \code{\link[spatstat]{envelope}}.
# @param verbose Logical flag indicating whether to print progress reports during the simulations.
# See \code{\link[spatstat]{envelope}}.
# @param MrkvickaEtal2017 Logical. If TRUE and type is "st" or "qdir", then the combined scaled MAD
# envelope presented in Mrkvička et al. (2017) is calculated. Otherwise, the two-step procedure
# described in \code{\link{global_envelope_test}} is used for combining the tests.
# @seealso \code{\link[spatstat]{envelope}} (that is used to perform simulations),
# \code{\link{global_envelope_test}}
#' @importFrom spatstat envelope
global_envelope_with_sims <- function(X, nsim, simfun = NULL, simfun.arg = NULL,
        testfuns = NULL, ...,
        type = "erl", alpha = 0.05, alternative = c("two.sided", "greater", "less"),
        r_min = NULL, r_max = NULL, take_residual = FALSE,
        ties = "erl", probs = c(0.025, 0.975),
        save.envelope = TRUE, savefuns = FALSE, savepatterns = FALSE,
        verbose = FALSE, MrkvickaEtal2017 = FALSE) {
  alt <- match.arg(alternative)
  nfuns <- length(testfuns)
  if(nfuns < 1) nfuns <- 1
  if(!is.null(r_min) & length(r_min) != nfuns) stop("r_min should give the minimum distances for each of the test functions.\n")
  if(!is.null(r_max) & length(r_max) != nfuns) stop("r_max should give the maximum distances for each of the test functions.\n")

  # Create simulations from the given model and calculate the test functions
  curve_sets_ls <- simpatterns_and_funcs_from_X(X, nsim=nsim, simfun=simfun, simfun.arg=simfun.arg,
                                                testfuns=testfuns, ...,
                                                savepatterns=savepatterns, verbose=verbose, calc_funcs=TRUE)
  # Crop curves (if r_min or r_max given)
  tmpfunc <- function(i) crop_curves(curve_sets_ls[[i]], r_min=r_min[i], r_max=r_max[i])
  curve_sets_ls <- lapply(1:nfuns, FUN=tmpfunc)

  # Make the test
  if(length(curve_sets_ls) > 1 & MrkvickaEtal2017 & type %in% c("st", "qdir")) {
    global_envtest <- combined_scaled_MAD_envelope(curve_sets_ls, type=type, alpha=alpha, probs=probs)
    stat <- attr(global_envtest$step2_test, "k")[1]
    pval <- attr(global_envtest$step2_test, "p")
    curve_set_combined <- attr(global_envtest, "rank_test_curve_set")
  }
  else {
    global_envtest <- global_envelope_test(curve_sets_ls, type=type, alpha=alpha,
                                           alternative=alt, ties=ties, probs=probs)
    switch(class(global_envtest)[1],
           global_envelope = {
             stat <- attr(global_envtest, "k")[1]
             pval <- attr(global_envtest, "p")
           },
           combined_global_envelope = {
             stat <- attr(global_envtest$step2_test, "k")[1]
             pval <- attr(global_envtest$step2_test, "p")
           })
  }

  res <- structure(list(statistic = as.numeric(stat), p.value = pval,
                  method = type, curve_set = curve_sets_ls), class = "global_envelope_with_sims")
  if(save.envelope) attr(res, "envelope_test") <- global_envtest
  if(savefuns) attr(res, "simfuns") <- curve_sets_ls # FIXME two times curve_set_ls saved
  if(savepatterns) attr(res, "simpatterns") <- attr(curve_sets_ls[[1]], "simpatterns")

  res
}

#' Adjusted global envelope tests
#'
#' Adjusted global envelope tests.
#'
#'
#' The specification of X is important here:
#'
#' 1) If simfun = NULL and fitfun = NULL (default), then \code{\link[spatstat]{envelope}}
#' is used for generating simulations under the null hypothesis and
#' \itemize{
#' \item If X is a point pattern, the null hypothesis is CSR.
#' \item If X is a fitted model, the null hypothesis is that model.
#' }
#'
#' 2) The user can provide the function for fitting the model (fitfun) and for simulating
#' from the fitted model (simfun). These functions should be coupled with each other such
#' that the object returned by 'fitfun' is directly accepted as the (single) argument in 'simfun'.
#' Further X should then be an \code{\link[spatstat]{ppp}} object and 'fitfun' should accept as
#' the argument an \code{\link[spatstat]{ppp}} object (X and further simulated point patterns).
#'
#'
#' A note: The structure of the code, which utilizes \code{\link[spatstat]{envelope}} though an
#' inner function, mimics the structure in the function \code{\link[spatstat]{dg.envelope}}
#' in the use of \code{\link[spatstat]{envelope}}. However, this function allows for more general
#' use as described above and for the different global envelope tests.
#'
#' For the rank envelope test, the test is the test described in Myllymäki et al. (2017) with the
#' adjustment of Baddeley et al. (2017).
#' For other test types, the test (also) uses the two-stage procedure of Dao and Genton (2014) with
#' the adjustment of Baddeley et al. (2017).
#'
#' See examples in \code{\link{saplings}}.
#'
#' @param X An object containing point pattern data. A point pattern (object of class "ppp")
#' or a fitted point process model (object of class "ppm" or "kppm"). See
#' \code{\link[spatstat]{envelope}}.
#' @param nsim The number of simulations to be generated in the primary test.
#' @param nsimsub Number of simulations in each basic test. There will be nsim repetitions of the
#' basic test, each involving nsimsub simulated realisations, so there will be a total of
#' nsim * (1 + nsimsub) simulations.
#' @param simfun A function for generating simulations from the null model. If given, this function
#' is called by replicate(n=nsim, simfun(simfun.param), simplify=FALSE) to make nsim simulations.
#' The function should return an \code{\link[spatstat]{ppp}} object as those are further passed to
#' \code{\link[spatstat]{envelope}}.
#' If the function is not provided, then \code{\link[spatstat]{envelope}} is used also for generating
#' the point patterns from the null hypothesis.
#' @param fitfun A function for estimating the parameters of the null model. If not given, then
#' \code{\link[spatstat]{envelope}} takes care of the parameter estimation as well (and X should be a fitted
#' model object). The function 'fitfun' should return the fitted model in the form that it can be directly
#' passed to 'simfun' as the argument 'simfun.arg'.
#' @param simfun.arg The parameter to be passed to simfun. The function simfun should handle
#' with the structure of simfun.param. FIXME should be here or not?
#' @param testfuns A list of lists of parameters to be passed to \code{\link[spatstat]{envelope}}.
#' A list of parameters should be provided for each test function that is to be used in the combined test.
#' @param ... Additional parameters to \code{\link[spatstat]{envelope}} in the case where only one test
#' function is used. In that case, this is an alternative to providing the parameter in the argument
#' testfuns. If \code{\link[spatstat]{envelope}} is also used to generate simulations under the null
#' hypothesis (if simfun not provided), then also recall to specify how to generate the simulations.
#' @inheritParams global_envelope_test
#' @param r_min The minimum radius to include in the test.
#' @param r_max The maximum radius to include in the test. Note: cannot be larger than r-values used
#' in calculating functions by \code{\link[spatstat]{envelope}}.
#' @param take_residual If (needed for visual reasons only) the theoretical or mean behaviour of the
#' test function is reduced from the test functions. If TRUE, then: If \code{\link[spatstat]{envelope}}
#' provides the theoretical value 'theo' for the simulated model, then this value is used. Otherwise,
#' the theoretical function is taken as the mean of the simulations.
# @param ties Ties method to be passed to \code{\link{global_envelope_test}}.
# Used to obtain a point estimate for the p-value for the "rank" test. Default to extreme rank
# length p-value.
#' @param save.cons.envelope Logical flag indicating whether to save the unadjusted envelope test results.
#' @param savefuns Logical flag indicating whether to save all the simulated function values.
#' See \code{\link[spatstat]{envelope}}.
#' @param savepatterns Logical flag indicating whether to save all the simulated point patterns.
#' See \code{\link[spatstat]{envelope}}.
#' @param verbose Logical flag indicating whether to print progress reports during the simulations.
#' See \code{\link[spatstat]{envelope}}.
#' @param MrkvickaEtal2017 Logical. If TRUE, type is "st" or "qdir" and several test functions are used,
#' then the combined scaled MAD envelope presented in Mrkvička et al. (2017) is calculated. Otherwise,
#' the two-step procedure described in \code{\link{global_envelope_test}} is used for combining the tests.
#' Default to FALSE.
#' @param mc.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously.
#' Must be at least one, and parallelization requires at least two cores. On a Windows computer mc.cores must be 1
#' (no parallelization). For details, see \code{\link[parallel]{mclapply}}, for which the argument is passed.
#'
#' @return An object of class \code{global_envelope} or \code{combined_global_envelope}.
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381-404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27(5): 1239-1255. DOI: 10.1007/s11222-016-9683-9
#'
#' Dao, N.A. and Genton, M. (2014). A Monte Carlo adjusted goodness-of-fit test for parametric models describing spatial point patterns. Journal of Graphical and Computational Statistics 23, 497-517.
#'
#' Baddeley, A., Hardegen, A., Lawrence, T., Milne, R. K., Nair, G. and Rakshit, S. (2017). On two-stage Monte Carlo tests of composite hypotheses. Computational Statistics and Data Analysis 114: 75-87. doi: http://dx.doi.org/10.1016/j.csda.2017.04.003
#' @seealso \code{\link{rank_envelope}}, \code{\link{qdir_envelope}}, \code{\link{st_envelope}},
#' \code{\link{plot.global_envelope}}, \code{\link{saplings}}
#' @export
#' @importFrom spatstat is.ppm
#' @importFrom spatstat is.kppm
#' @importFrom spatstat is.lppm
#' @importFrom spatstat is.slrm
#' @importFrom spatstat is.ppp
#' @importFrom spatstat update.ppm
#' @importFrom spatstat update.kppm
#' @importFrom spatstat update.lppm
#' @importFrom spatstat update.slrm
#' @importFrom parallel mclapply
#' @importFrom stats quantile
dg.global_envelope_test <- function(X, nsim = 499, nsimsub = nsim,
        simfun=NULL, fitfun=NULL, simfun.arg=NULL, testfuns=NULL, ..., type = "erl",
        alpha = 0.05, alternative = c("two.sided","less", "greater"),
        r_min=NULL, r_max=NULL, take_residual=FALSE,
        save.cons.envelope = savefuns || savepatterns, savefuns = FALSE, savepatterns = FALSE,
        verbose=TRUE, MrkvickaEtal2017 = FALSE, mc.cores=1L) {
    Xismodel <- spatstat::is.ppm(X) || spatstat::is.kppm(X) || spatstat::is.lppm(X) || spatstat::is.slrm(X)
    if(is.null(simfun) != is.null(fitfun)) stop("Either both \'simfun\' and \'fitfun\' should be provided or neither of them.\n")
    if((is.null(simfun) | is.null(fitfun)) & !Xismodel)
        warning("\'simfun\' and/or \'fitfun\' not provided and \'X\' is not a fitted model object.\n",
                "As \'envelope\' is used for simulations and model fitting, \n",
                "\'X\' should be a fitted model object.")
    if(!is.null(fitfun) & !spatstat::is.ppp(X)) stop("Model to be fitted by fitfun(X) and simulations to be generated by simfun,\n but X is not an ppp object.\n")
    nfuns <- length(testfuns)
    if(nfuns < 1) nfuns <- 1
    alt <- match.arg(alternative)
    simfun.arg <- NULL
    if(!is.null(fitfun)) simfun.arg <- fitfun(X) # fitted model parameters to be passed to simfun
    if(verbose) cat("Applying test to original data...\n")
    tX <- global_envelope_with_sims(X, nsim=nsim, simfun=simfun, simfun.arg=simfun.arg,
            testfuns=testfuns, ..., type = type, alpha = alpha, alternative = alt,
            r_min=r_min, r_max=r_max, take_residual=take_residual,
            ties='midrank', # Note: ties argument do not matter here (p-values not used for the rank test), using 'midrank' as a faster alternative.
            save.envelope = save.cons.envelope, savefuns = TRUE, savepatterns = savepatterns,
            verbose = verbose, MrkvickaEtal2017=MrkvickaEtal2017)
    if(verbose) cat("Done.\n")
    # Create another set of simulated patterns to be used to estimate the second-state p-value
    # following Baddeley et al. (2017).
    # FIXME: unnecessary calculations of test functions in simpatterns_and_funcs_from_X
    simpatterns <- attr(simpatterns_and_funcs_from_X(X, nsim=nsim, simfun=simfun, simfun.arg=simfun.arg, testfuns=testfuns, ...,
                                                     savepatterns=TRUE, verbose=verbose,
                                                     calc_funcs=FALSE)[[1]], "simpatterns")

    if(verbose) cat(paste("Running tests on", nsim, "simulated patterns... \n"))
    # For each of the simulated patterns in 'simpatterns', perform the test and calculate
    # the extreme rank (or deviation) measure and p-value
    loopfun <- function(rep) {
      Xsim <- simpatterns[[rep]]
      if(!is.null(fitfun)) simfun.arg <- fitfun(Xsim) # a fitted model to be passed to simfun
      if(Xismodel)
          switch(class(X)[1],
                 ppm = {
                     Xsim <- spatstat::update.ppm(X, Xsim)
                 },
                 kppm = {
                     Xsim <- spatstat::update.kppm(X, Xsim)
                 },
                 lppm = {
                     Xsim <- spatstat::update.lppm(X, Xsim)
                 },
                 slrm = {
                     Xsim <- spatstat::update.slrm(X, Xsim)
                 })
      tXsim <- global_envelope_with_sims(Xsim, nsim=nsimsub, simfun=simfun, simfun.arg=simfun.arg,
              testfuns=testfuns, ..., type = type, alpha = alpha, alternative = alt,
              r_min=r_min, r_max=r_max, take_residual=take_residual,
              ties='midrank', # Note: the ties method does not matter here
              save.envelope = FALSE, savefuns = FALSE, savepatterns = FALSE,
              verbose = verbose, MrkvickaEtal2017=MrkvickaEtal2017)
      list(stat = tXsim$statistic, pval = tXsim$p.value)
    }
    mclapply_res <- parallel::mclapply(X=1:nsim, FUN=loopfun, mc.cores=mc.cores)
    stats <- sapply(mclapply_res, function(x) x$stat)
    pvals <- sapply(mclapply_res, function(x) x$pval)

    # The adjusted test
    if(MrkvickaEtal2017 & type %in% c("st", "qdir") & length(nfuns) > 1) { # Combined tests as in Mrkvicka et al. (2017)
      res <- attr(tX, "envelope_test")
      #-- The rank test at the second level
      # Calculate the critical rank / alpha
      kalpha_star <- stats::quantile(stats, probs=alpha, type=1)
      data_and_sim_curves <- data_and_sim_curves(attr(res, "level2_curve_set")) # the second step "curves"
      nr <- ncol(data_and_sim_curves)
      Nfunc <- nrow(data_and_sim_curves)
      LB <- array(0, nr)
      UB <- array(0, nr)
      for(i in 1:nr){
        Hod <- sort(data_and_sim_curves[,i])
        LB[i]<- Hod[kalpha_star]
        UB[i]<- Hod[Nfunc-kalpha_star+1]
      }
      # -> kalpha_star, LB, UB of the (second level) rank test
      # Update res object with adjusted values
      res$lo <- LB
      res$hi <- UB
      attr(res, "method") <- "Adjusted global envelope test" # Change method name
      attr(res, "k_alpha_star") <- kalpha_star # Add kalpha_star
      # Re-calculate the new qdir/st envelopes
      envchars <- combined_scaled_MAD_bounding_curves_chars(attr(tX, "simfuns"), type = type)
      central_curves_ls <- lapply(attr(tX, "simfuns"), function(x) get_T_0(x))
      bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=UB,
                                                             lower_f=envchars$lower_f, upper_f=envchars$upper_f)
      # Update the first level envelopes for plotting
      for(i in 1:length(attr(res, "global_envelope_ls"))) {
        attr(res, "global_envelope_ls")[[i]]$lo <- bounding_curves$lower_ls[[i]]
        attr(res, "global_envelope_ls")[[i]]$hi <- bounding_curves$upper_ls[[i]]
      }
    }
    else { # Otherwise, the ERL test is used at the second level of a combined test
      if(type == "rank" & nfuns == 1) {
        # Calculate the critical rank (instead of alpha) and the adjusted envelope following Myllymäki et al. (2017)
        kalpha_star <- stats::quantile(stats, probs=alpha, type=1)
        data_and_sim_curves <- data_and_sim_curves(tX$curve_set[[1]]) # all the functions
        nr <- length(tX$curve_set[[1]]$r)
        Nfunc <- nrow(data_and_sim_curves)
        LB <- array(0, nr)
        UB <- array(0, nr)
        for(i in 1:nr){
          Hod <- sort(data_and_sim_curves[,i])
          LB[i]<- Hod[kalpha_star]
          UB[i]<- Hod[Nfunc-kalpha_star+1]
        }
        # For getting automatically an global_envelope object, first call central_region
        res <- central_region(tX$curve_set, type=type, coverage=1-alpha, alternative=alt, central="mean")
        # Update res object with adjusted values
        res$lo <- LB
        res$hi <- UB
        attr(res, "method") <- "Adjusted global envelope test" # Change method name
        attr(res, "k_alpha_star") <- kalpha_star # Add kalpha_star
      }
      else {
        alpha_star <- stats::quantile(pvals, probs=alpha, type=1)
        res <- global_envelope_test(tX$curve_set, type=type, alpha=alpha_star, alternative=alt)
        # Save additional arguments
        attr(res, "method") <- "Adjusted global envelope test" # Change method name
        attr(res, "alpha") <- alpha
        attr(res, "alpha_star") <- alpha_star
      }
    }

    # Changes
    attr(res, "call") <- match.call() # Update "call" attribute
    # Additions
    if(savepatterns) attr(res, "simpatterns") <- simpatterns
    if(savefuns) attr(res, "simfuns") <- attr(tX, "simfuns")
    if(save.cons.envelope) attr(res, "unadjusted_envelope_test") <- attr(tX, "envelope_test")
    attr(res, "simulated_p") <- pvals
    attr(res, "simulated_k") <- stats
    # Return
    res
}
