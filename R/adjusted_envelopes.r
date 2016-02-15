#' A global envelope test
#'
#' A global envelope test including simulations from a point process model.
#'
#'
#' Details of the tests are given in \code{\link{rank_envelope}} (rank envelope test),
#' \code{\link{qdir_envelope}} (directional quantile envelope test) and \code{\link{st_envelope}}
#' (studentized envelope test).
#'
#' The specification of X is important here:
#' if simfun is not provided, the function \code{\link[spatstat]{envelope}} is used to generate
#' simulations under the null hypothesis and to calculate the test functions (specified in the
#' arguments ...) and then
#' \itemize{
#' \item If X is a point pattern, the null hypothesis is CSR.
#' \item If X is a fitted model, the null hypothesis is that model.
#' }
#' If simfun is provided, then the null model is the one simulated by this given function,
#' and X is expected to be a point pattern of \code{\link[spatstat]{ppp}} object, for which data
#' pattern and simulations \code{\link[spatstat]{envelope}} calculates the test statistics.
#'
#' @param X An object containing point pattern data. A point pattern (object of class "ppp")
#' or a fitted point process model (object of class "ppm" or "kppm"). See
#' \code{\link[spatstat]{envelope}}.
#' @param nsim The number of simulations.
#' @param simfun A function for generating simulations from the null model. If given, this function
#' is called by replicate(n=nsim, simfun(simfun.param), simplify=FALSE) to make nsim simulations.
#' The function should return an \code{\link[spatstat]{ppp}} object as those are further passed to
#' \code{\link[spatstat]{envelope}}.
#' If the function is not provided, then \code{\link[spatstat]{envelope}} is used also for generating
#' the point patterns from the null hypothesis.
#' @param simfun.arg The parameter to be passed to simfun. The function simfun should handle
#' with the structure of simfun.param.
#' @param ... Additional parameters passed to \code{\link[spatstat]{envelope}}.
#' For example, the test function in the argument 'fun' and further specifications regarding that.
#' If \code{\link[spatstat]{envelope}} is also used to generate simulations under the null hypothesis
#' (if simfun not provided), then also recall to specify how to generate the simulations.
#' @param test Either "rank" for the \code{\link{rank_envelope}} test, "qdir" for the
#' \code{\link{qdir_envelope}} test or "st" for the \code{\link{st_envelope}} test.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of
#' the following: "two.sided" (default), "less" or "greater" for "rank". Relevant only for the
#' rank test (otherwise ignored).
#' @param r_min The minimum radius to include in the test.
#' @param r_max The maximum radius to include in the test. Note: cannot be larger than r-values used
#' in calculating functions by \code{\link[spatstat]{envelope}}.
#' @param take_residual If (needed for visual reasons only) the theoretical or mean behaviour of the
#' test function is reduced from the test functions. If TRUE, then: If \code{\link[spatstat]{envelope}}
#' provides the theoretical value 'theo' for the simulated model, then this value is used. Otherwise,
#' the theoretical function is taken as the mean of the simulations.
#' @param lexo Logical, whether or not to use lexical ordering for calculation of the p-value.
#' in the rank envelope test. See \code{\link{rank_envelope}}.
#' @param ties Ties method to be passed to \code{\link{rank_envelope}} (and then to
#' \code{\link{estimate_p_value}}. Used to obtain a point estimate for the p-value, if lexo=FALSE.
#' The default point estimate is the mid-rank p-value.
#' @param save.envelope Logical flag indicating whether to save all envelope test results.
#' @param savefuns Logical flag indicating whether to save all the simulated function values.
#' See \code{\link[spatstat]{envelope}}.
#' @param savepatterns Logical flag indicating whether to save all the simulated point patterns.
#' See \code{\link[spatstat]{envelope}}.
#' @param verbose Logical flag indicating whether to print progress reports during the simulations.
#' See \code{\link[spatstat]{envelope}}.
#' @seealso \code{\link[spatstat]{envelope}} (that is used to perform simulations),
#' \code{\link{rank_envelope}}, \code{\link{qdir_envelope}}, \code{\link{st_envelope}}
global_envelope_with_sims <- function(X, nsim, simfun=NULL, simfun.arg=NULL, ..., test = c("rank", "qdir", "st"),
        alpha = 0.05, alternative = c("two.sided", "less", "greater"),
        r_min=NULL, r_max=NULL, take_residual=FALSE,
        lexo = TRUE, ties=NULL,
        save.envelope = TRUE, savefuns = FALSE, savepatterns = FALSE,
        verbose = FALSE) {
    test <- match.arg(test)
    alt <- match.arg(alternative)
    if(!is.null(simfun)) {
        # Create simulations by the given function
        simpatterns <- replicate(n=nsim, simfun(simfun.arg), simplify=FALSE)
        # Calculate the test functions
        X <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ..., savefuns = TRUE, savepatterns = savepatterns, verbose=verbose)
    }
    else {
        # Create simulations from the given model and calculate the test functions
        X <- spatstat::envelope(X, nsim=nsim, ..., savefuns = TRUE, savepatterns = savepatterns, verbose=verbose)
    }
    # Crop curves (if r_min or r_max given)
    curve_set <- crop_curves(X, r_min = r_min, r_max = r_max)
    # Make the test
    switch(test,
           rank = {
               global_envtest <- rank_envelope(curve_set, alpha=alpha, savedevs=TRUE,
                       alternative=alt, lexo=lexo, ties=ties)
               stat <- global_envtest$k[1]
           },
           qdir = {
               global_envtest <- qdir_envelope(curve_set, alpha=alpha, savedevs=TRUE)
               stat <- global_envtest$u[1]
           },
           st = {
               global_envtest <- st_envelope(curve_set, alpha=alpha, savedevs=TRUE)
               stat <- global_envtest$u[1]
           })

    res <- structure(list(statistic = as.numeric(stat), p.value = global_envtest$p,
                          method = test, curve_set = curve_set), class = "global_envelope_with_sims")
    if(save.envelope) {
        attr(res, "envelope_test") <- global_envtest
    }
    if(savefuns) attr(res, "simfuns") <- attr(X, "simfuns")
    if(savepatterns) attr(res, "simpatterns") <- attr(X, "simpatterns")

    res
}


#' A combined global envelope test
#'
#' A combined global envelope test including simulations from a point process model.
#'
#'
#' Details of the tests are given in \code{\link{rank_envelope}} (rank envelope test),
#' \code{\link{combined_scaled_MAD_envelope}} (directional quantile and studentized envelope tests).
#'
#' The specification of X is important here:
#' if simfun is not provided, the function \code{\link[spatstat]{envelope}} is used to generate
#' simulations under the null hypothesis and to calculate the test functions (specified in the
#' arguments ...) and then
#' \itemize{
#' \item If X is a point pattern, the null hypothesis is CSR.
#' \item If X is a fitted model, the null hypothesis is that model.
#' }
#' If simfun is provided, then the null model is the one simulated by this given function,
#' and X is expected to be a point pattern of \code{\link[spatstat]{ppp}} object, for which data
#' pattern and simulations \code{\link[spatstat]{envelope}} calculates the test statistics.
#'
#' If savefuns is TRUE, all the simulated functions are saved in an attribute "simfuns" as a list
#' of curve sets for each test function.
#' @inheritParams global_envelope_with_sims
#' @seealso \code{\link[spatstat]{envelope}} (that is used to perform simulations),
#' \code{\link{rank_envelope}}, \code{\link{qdir_envelope}}, \code{\link{st_envelope}}
combined_global_envelope_with_sims <- function(X, nsim, simfun=NULL, simfun.arg=NULL, ...,
        testfuns = NULL,
        test = c("qdir", "st", "rank"),
        alpha = 0.05, alternative = c("two.sided", "less", "greater"),
        r_min=NULL, r_max=NULL, take_residual=FALSE,
        lexo = TRUE, ties=NULL,
        save.envelope = TRUE, savefuns = FALSE, savepatterns = FALSE,
        verbose = FALSE) {
    test <- match.arg(test)
    alt <- match.arg(alternative)

    nfuns <- length(testfuns)
    if(length(r_min) != nfuns) stop("r_min should give the minimum distances for each of the test functions.\n")
    if(length(r_max) != nfuns) stop("r_max should give the maximum distances for each of the test functions.\n")
    testfuns_argnames_ls <- lapply(testfuns, function(x) names(x))
    if(!is.null(simfun)) {
        # Create simulations by the given function
        simpatterns <- replicate(n=nsim, simfun(simfun.arg), simplify=FALSE)
        # Calculate the test functions
        for(i in 1:nfuns) {
            tmpstring <- paste("X.testfunc", i, " <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ", sep="")
            for(j in 1:length(testfuns_argnames_ls[[i]]))
                tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], " = testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
            tmpstring <- paste(tmpstring, "savefuns = TRUE, savepatterns = FALSE, verbose=verbose, ...)", sep="")
            eval(parse(text = tmpstring))
        }
    }
    else {
        # Create simulations from the given model and calculate the test functions
        i = 1
        tmpstring <- paste("X.testfunc", i, " <- spatstat::envelope(X, nsim=nsim, ", sep="")
        for(j in 1:length(testfuns_argnames_ls[[i]]))
            tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], " = testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
        tmpstring <- paste(tmpstring, "savefuns = TRUE, savepatterns = TRUE, verbose=verbose, ...)", sep="")
        eval(parse(text = tmpstring))
        simpatterns <- attr(X.testfunc1, "simpatterns")
        for(i in 2:nfuns) {
            tmpstring <- paste("X.testfunc", i, " <- spatstat::envelope(X, nsim=nsim, simulate=simpatterns, ", sep="")
            for(j in 1:length(testfuns_argnames_ls[[i]]))
                tmpstring <- paste(tmpstring, paste(testfuns_argnames_ls[[i]][j], " = testfuns[[", i, "]]", "[[", j, "]], ", sep=""), sep="")
            tmpstring <- paste(tmpstring, "savefuns = TRUE, savepatterns = FALSE, verbose=verbose, ...)", sep="")
            eval(parse(text = tmpstring))
        }
    }
    curve_sets_ls <- NULL
    for(i in 1:nfuns) curve_sets_ls[[i]] <- get(paste("X.testfunc", i, sep=""))
    # Crop curves (if r_min or r_max given)
    tmpfunc <- function(i) crop_curves(curve_sets_ls[[i]], r_min = r_min[i], r_max = r_max[i])
    curve_sets_ls <- lapply(1:nfuns, FUN = tmpfunc)

    # Make the test
    switch(test,
            rank = {
                # Create a combined curve_set
                curve_set_combined <- combine_curve_sets(curve_sets_ls)
                global_envtest <- rank_envelope(curve_set_combined, alpha=alpha, savedevs=TRUE,
                        alternative=alt, lexo=lexo, ties=ties)
                stat <- global_envtest$k[1]
                pval <- global_envtest$p
            },
            qdir = {
                global_envtest <- combined_scaled_MAD_envelope(curve_sets_ls, test="qdir", alpha=alpha)
                stat <- global_envtest$rank_test$k[1]
                pval <- global_envtest$rank_test$p
                curve_set_combined <- attr(global_envtest, "rank_test_curve_set")
            },
            st = {
                global_envtest <- combined_scaled_MAD_envelope(curve_sets_ls, test="st", alpha=alpha)
                stat <- global_envtest$rank_test$k[1]
                pval <- global_envtest$rank_test$p
                curve_set_combined <- attr(global_envtest, "rank_test_curve_set")
            })

    res <- structure(list(statistic = as.numeric(stat), p.value = pval,
                    method = test, curve_set = curve_set_combined), class = "global_envelope_with_sims")
    if(save.envelope) attr(res, "envelope_test") <- global_envtest
    if(savefuns) attr(res, "simfuns") <- curve_sets_ls
    if(savepatterns) attr(res, "simpatterns") <- simpatterns

    res
}

#' Adjusted global envelope tests
#'
#' Adjusted global rank envelope test, studentized envelope test and directional quantile envelope test.
#'
#'
#' The specification of X is important here:
#' \itemize{
#' \item If X is a point pattern, the null hypothesis is CSR.
#' \item If X is a fitted model, the null hypothesis is that model.
#' }
#'
#' The structure of the code, which utilizes \code{\link[spatstat]{envelope}} though the function
#' \code{\link{global_envelope_with_sims}}, mimics the structure in the function
#' \code{\link[spatstat]{dg.envelope}} in the use of \code{\link[spatstat]{envelope}}.
#'
#' @inheritParams global_envelope_with_sims
#' @param nsim The number of simulations to be generated in the primary test.
#' @param nsimsub Number of simulations in each basic test. There will be nsim repetitions of the
#' basic test, each involving nsimsub simulated realisations, so there will be a total of
#' nsim * (1 + nsimsub) simulations.
#' @param fitfun A function for estimating the parameters of the null model. If not given, then
#' \code{\link[spatstat]{envelope}} takes care of the parameter estimation as well (and X should be a fitted
#' model object). The function 'fitfun' should return the fitted model in the form that it can be directly
#' passed to 'simfun' as the argument 'simfun.arg'.
#' @param save.cons.envelope Logical flag indicating whether to save the unadjusted envelope test results.
#' @param mc.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously.
#' Must be at least one, and parallelization requires at least two cores. On a Windows computer mc.cores must be 1
#' (no parallelization). For details, see \code{\link[parallel]{mclapply}}, for which the argument is passed.
#'
#' @return An object of class adjusted_envelope_test.
#' @references
#' Dao, N.A. and Genton, M. (2014). A Monte Carlo adjusted goodness-of-fit test for parametric models describing spatial point patterns. Journal of Graphical and Computational Statistics 23, 497-517.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2015). Global envelope tests for spatial point patterns. arXiv:1307.0239v4 [stat.ME]
#'
#' @seealso \code{\link{rank_envelope}}, \code{\link{qdir_envelope}}, \code{\link{st_envelope}},
#' \code{\link{plot.adjusted_envelope_test}}
#' @export
dg.global_envelope <- function(X, nsim = 499, nsimsub = nsim,
        simfun=NULL, fitfun=NULL, ..., test = c("rank", "qdir", "st"),
        alpha = 0.05, alternative = c("two.sided","less", "greater"),
        r_min=NULL, r_max=NULL, take_residual=FALSE,
        #rank_count_test_p_values = FALSE, lexo = TRUE, ties=NULL,
        save.cons.envelope = savefuns || savepatterns, savefuns = FALSE, savepatterns = FALSE,
        verbose=TRUE, mc.cores=1L) {
    Xismodel <- spatstat::is.ppm(X) || spatstat::is.kppm(X) || spatstat::is.lppm(X) || spatstat::is.slrm(X)
    if(is.null(simfun) != is.null(fitfun)) stop("Either both \'simfun\' and \'fitfun\' should be provided or neither of them.\n")
    if((is.null(simfun) | is.null(fitfun)) & !Xismodel)
        warning("\'simfun\' and/or \'fitfun\' not provided and \'X\' is not a fitted model object.\n",
                "As \'envelope\' is used for simulations and model fitting, \n",
                "\'X\' should be a fitted model object.")
    test <- match.arg(test)
    alt <- match.arg(alternative)
    simfun.arg <- NULL
    if(!is.null(fitfun)) simfun.arg <- fitfun(X) # a fitted model parameters to be passed to simfun
    if(verbose) cat("Applying test to original data...\n")
    tX <- global_envelope_with_sims(X, nsim=nsim, simfun=simfun, simfun.arg=simfun.arg, ...,
            test = test, alpha = alpha, alternative = alt,
            r_min=r_min, r_max=r_max, take_residual=take_residual,
            lexo = FALSE, ties='midrank',
            save.envelope = save.cons.envelope, savefuns = savefuns, savepatterns = TRUE,
            verbose = verbose)
    if(verbose) cat("Done.\n")
    simpatterns <- attr(tX, "simpatterns")

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
        tXsim <- global_envelope_with_sims(Xsim, nsim=nsimsub, simfun=simfun, simfun.arg=simfun.arg, ...,
                test = test, alpha = alpha, alternative = alt,
                r_min=r_min, r_max=r_max, take_residual=take_residual,
                lexo = FALSE, ties='midrank', # Note: the ties method does not matter here; p-values not used for the rank test.
                save.envelope = FALSE, savefuns = FALSE, savepatterns = FALSE,
                verbose = verbose)
        list(stat = tXsim$statistic, pval = tXsim$p.value)
    }
    mclapply_res <- parallel::mclapply(X=1:nsim, FUN=loopfun, mc.cores=mc.cores)
    stats <- sapply(mclapply_res, function(x) x$stat)
    pvals <- sapply(mclapply_res, function(x) x$pval)

    # Calculate the critical rank / alpha
    switch(test,
           rank = {
               kalpha_star <- quantile(stats, probs=alpha, type=1)
               #alpha_star <- quantile(pvals, probs=alpha, type=1)

               data_curve <- tX$curve_set[['obs']]
               sim_curves <- t(tX$curve_set[['sim_m']])
               data_and_sim_curves <- rbind(data_curve, sim_curves)
               T_0 <- get_T_0(tX$curve_set)

               nr <- length(tX$curve_set$r)
               LB <- array(0, nr);
               UB <- array(0, nr);
               for(i in 1:nr){
                   Hod <- sort(data_and_sim_curves[,i])
                   LB[i]<- Hod[kalpha_star]
                   UB[i]<- Hod[nsim+1-kalpha_star+1]
               }

               adjenv <- list(r=tX$curve_set[['r']], method="Adjusted rank envelope test",
                       alternative = alt,
                       p=NULL, p_interval=NULL, 
                       k_alpha=kalpha_star, k=stats,
                       central_curve=T_0, data_curve=data_curve, lower=LB, upper=UB,
                       call=match.call())
               test_name <- "Adjusted rank envelope test"
           },
           qdir = {
               alpha_star <- quantile(pvals, probs=alpha, type=1)
               adjenv <- qdir_envelope(tX$curve_set, alpha=alpha_star)
               adjenv$alpha_star <- alpha_star
               adjenv$p_values <- pvals
               adjenv$alpha <- alpha
               test_name <- "Adjusted directional quantile envelope test"
           },
           st = {
               alpha_star <- quantile(pvals, probs=alpha, type=1)
               adjenv <- st_envelope(tX$curve_set, alpha=alpha_star)
               adjenv$alpha_star <- alpha_star
               adjenv$p_values <- pvals
               adjenv$alpha <- alpha
               test_name <- "Adjusted studentized envelope test"
           })

    res <- structure(list(adj_envelope_test = adjenv,
                    method = test_name), class = "adjusted_envelope_test")
    if(savepatterns) attr(res, "simpatterns") <- simpatterns
    if(savefuns) attr(res, "simfuns") <- attr(tX, "simfuns")
    if(save.cons.envelope) attr(res, "unadjusted_envelope_test") <- attr(tX, "envelope_test")
    attr(res, "alternative") <- alt
    attr(res, "call") <- match.call()
    res
}


#' Adjusted combined global envelope tests
#'
#' Adjusted combined global rank envelope test, studentized envelope test or directional quantile envelope test.
#'
#'
#' The specification of X is important here:
#' \itemize{
#' \item If X is a point pattern, the null hypothesis is CSR.
#' \item If X is a fitted model, the null hypothesis is that model.
#' }
#'
#' The structure of the code, which utilizes \code{\link[spatstat]{envelope}} though the function
#' \code{\link{global_envelope_with_sims}}, mimics the structure in the function
#' \code{\link[spatstat]{dg.envelope}} in the use of \code{\link[spatstat]{envelope}}.
#'
#' @inheritParams combined_global_envelope_with_sims
#' @inheritParams dg.global_envelope
#'
#' @return An object of class adjusted_envelope_test.
#' @references
#' Mrkvicka, T., Myllymäki, M. and Hahn, U. (2015). Multiple Monte Carlo testing with applications in spatial point processes. arXiv:1506.01646 [stat.ME]
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2015). Global envelope tests for spatial point patterns. arXiv:1307.0239v4 [stat.ME]
#'
#' Dao, N.A. and Genton, M. (2014). A Monte Carlo adjusted goodness-of-fit test for parametric models describing spatial point patterns. Journal of Graphical and Computational Statistics 23, 497-517.
#'
#' @seealso \code{\link{rank_envelope}}, \code{\link{qdir_envelope}}, \code{\link{st_envelope}},
#' \code{\link{plot.adjusted_envelope_test}}
#' @export
dg.combined_global_envelope <- function(X, nsim = 499, nsimsub = nsim,
        simfun=NULL, fitfun=NULL, ...,
        testfuns = NULL, test = c("qdir", "st", "rank"),
        alpha = 0.05, alternative = c("two.sided","less", "greater"),
        r_min=NULL, r_max=NULL, take_residual=FALSE,
        #rank_count_test_p_values = FALSE, lexo = TRUE, ties=NULL,
        save.cons.envelope = savefuns || savepatterns, savefuns = FALSE, savepatterns = FALSE,
        verbose=TRUE, mc.cores=1L) {
    Xismodel <- spatstat::is.ppm(X) || spatstat::is.kppm(X) || spatstat::is.lppm(X) || spatstat::is.slrm(X)
    if(is.null(simfun) != is.null(fitfun)) stop("Either both \'simfun\' and \'fitfun\' should be provided or neither of them.\n")
    if((is.null(simfun) | is.null(fitfun)) & !Xismodel)
        warning("\'simfun\' and/or \'fitfun\' not provided and \'X\' is not a fitted model object.\n",
                "As \'envelope\' is used for simulations and model fitting, \n",
                "\'X\' should be a fitted model object.")
    test <- match.arg(test)
    alt <- match.arg(alternative)
    simfun.arg <- NULL
    if(!is.null(fitfun)) simfun.arg <- fitfun(X) # a fitted model parameters to be passed to simfun
    if(verbose) cat("Applying test to original data...\n")
    tX <- combined_global_envelope_with_sims(X, nsim=nsim, simfun=simfun, simfun.arg=simfun.arg, ...,
            testfuns = testfuns, test = test, alpha = alpha, alternative = alt,
            r_min=r_min, r_max=r_max, take_residual=take_residual,
            lexo = FALSE, ties='midrank',
            save.envelope = save.cons.envelope, savefuns = savefuns, savepatterns = TRUE,
            verbose = verbose)
    if(verbose) cat("Done.\n")
    simpatterns <- attr(tX, "simpatterns")

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
        tXsim <- combined_global_envelope_with_sims(Xsim, nsim=nsimsub, simfun=simfun, simfun.arg=simfun.arg, ...,
                testfuns = testfuns, test = test, alpha = alpha, alternative = alt,
                r_min=r_min, r_max=r_max, take_residual=take_residual,
                lexo = FALSE, ties='midrank', # Note: the ties method does not matter here; p-values not used for the rank test.
                save.envelope = FALSE, savefuns = FALSE, savepatterns = FALSE,
                verbose = verbose)
        list(stat = tXsim$statistic, pval = tXsim$p.value)
    }
    mclapply_res <- parallel::mclapply(X=1:nsim, FUN=loopfun, mc.cores=mc.cores)
    stats <- sapply(mclapply_res, function(x) x$stat)
    pvals <- sapply(mclapply_res, function(x) x$pval)

    #-- The rank test
    # Calculate the critical rank / alpha
    kalpha_star <- quantile(stats, probs=alpha, type=1)
    #alpha_star <- quantile(pvals, probs=alpha, type=1)

    data_curve <- tX$curve_set[['obs']]
    sim_curves <- t(tX$curve_set[['sim_m']])
    data_and_sim_curves <- rbind(data_curve, sim_curves)
    T_0 <- get_T_0(tX$curve_set)

    nr <- length(tX$curve_set$r)
    LB <- array(0, nr);
    UB <- array(0, nr);
    for(i in 1:nr){
        Hod <- sort(data_and_sim_curves[,i])
        LB[i]<- Hod[kalpha_star]
        UB[i]<- Hod[nsim+1-kalpha_star+1]
    }
    # -> kalpha_stat, LB, UB of the rank test
    adjenv <- list(r=tX$curve_set[['r']], method="Adjusted rank envelope test",
            alternative = alt,
            p=NULL, p_interval=NULL,
            k_alpha=kalpha_star, k=stats,
            central_curve=T_0, data_curve=data_curve, lower=LB, upper=UB,
            call=match.call())
    test_name <- "Adjusted rank envelope test"

    # In the case of the combined scaled MAD envelope tests, we also need to calculate the new qdir/st envelopes
    res_env <- NULL
    if(test != "rank") {
        envchars <- combined_scaled_MAD_bounding_curves_chars(attr(tX, "simfuns"), test = test)
        central_curves_ls <- lapply(attr(tX, "simfuns"), function(x) get_T_0(x))
        bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=adjenv$upper,
                                                               lower_f=envchars$lower_f, upper_f=envchars$upper_f)
        # Create a combined envelope object for plotting
        res_env <- structure(list(r = do.call(c, lapply(attr(tX, "simfuns"), FUN = function(x) x$r), quote=FALSE),
                                  method = "Adjusted combined scaled MAD envelope test",
                                  alternative = "two.sided",
                                  p = adjenv$p,
                                  p_interval = adjenv$p_interval,
                                  central_curve = as.vector(do.call(c, central_curves_ls, quote=FALSE)),
                                  data_curve = do.call(c, lapply(attr(tX, "simfuns"), FUN = function(x) x[['obs']]), quote=FALSE),
                                  lower = do.call(c, bounding_curves$lower_ls, quote=FALSE),
                                  upper = do.call(c, bounding_curves$upper_ls, quote=FALSE)),
                             class = c("combined_scaled_MAD_envelope", "envelope_test"))
    }

    res <- structure(list(adj_envelope_test = adjenv,
                          method = test_name), class = "adjusted_envelope_test")
    if(savepatterns) attr(res, "simpatterns") <- simpatterns
    if(savefuns) attr(res, "simfuns") <- attr(tX, "simfuns")
    if(save.cons.envelope) attr(res, "unadjusted_envelope_test") <- attr(tX, "envelope_test")
    if(!is.null(res_env)) attr(res, "adjusted_scaled_MAD_envelope") <- res_env
    attr(res, "alternative") <- alt
    attr(res, "call") <- match.call()
    res
}


#' Print method for the class 'envelope_test'
#' @usage \method{print}{adjusted_envelope_test}(x, ...)
#'
#' @param x an 'adjusted_envelope_test' object
#' @param ... Ignored.
#'
#' @method print adjusted_envelope_test
#' @export
print.adjusted_envelope_test <- function (x, ...) {
    cat("Plot the object instead.\n")
}

#' Plot method for the class 'adjusted_envelope_test'
#' @usage \method{plot}{adjusted_envelope_test}(x, main,
#' plot_unadjusted=!is.null(attr(x, "unadjusted_envelope_test")), ...)
#'
#' @param x an 'adjusted_envelope_test' object
#' @param main See \code{\link{plot.default}}. Default is x$method.
#' @param plot_unadjusted Logical whether or not to plot also the unadjusted envelope.
#' Only available if these have been saved in 'x'.
#' @param ... Additional parameters to be passed to \code{\link{plot.envelope_test}}, if plot_unadjusted
#' is FALSE, or to \code{\link{two_envelopes_ggplot}}, if plot_unadjusted is TRUE.
#'
#' @method plot adjusted_envelope_test
#' @export
#' @seealso \code{\link{plot.envelope_test}}, \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}.
plot.adjusted_envelope_test <- function (x, main,
        plot_unadjusted=!is.null(attr(x, "unadjusted_envelope_test")), ...) {
    if(missing(main)) main <- x$method
    if(!plot_unadjusted) p <- plot.envelope_test(x$adj_envelope_test, main=main, ...)
    else {
        two_envelopes_ggplot(env1 = x$adj_envelope_test, env2 = attr(x, "unadjusted_envelope_test"))
    }
}
