#' The rank envelope test
#'
#' The rank envelope test, p-value and global envelope
#'
#'
#' The rank envelope test is a completely non-parametric test, which provides
#' the 100(1-alpha)\% global envelope for the chosen test function T(r) on
#' the chosen interval of distances and a p-value interval given by the most
#' liberal and the most conservative p-value estimate.
#'
#' Given a \code{curve_set} (or an \code{\link[spatstat]{envelope}}) object,
#' which contains both the data curve T_1(r) and the simulated curves T_2(r),...T_(s+1)(r),
#' the test is carried out as follows.
#'
#' For each curve in the curve_set, both the data curve and the simulations,
#' the global rank measure R is determined. If savedevs = TRUE, then the
#' global rank values R_1, R_2, ..., R_(s+1) are returned in the component 'k',
#' where k[1] is the value for the data.
#'
#' Based on R_i, i=1, ..., s+1, the p-interval is calculated. This interval is
#' by default plotted for the object returned by the rank_envelope function.
#' Also a single p-value is calculated and returned in component 'p'. By default
#' this p-value is the mid-rank p-value, but another option can be used by specifying
#' \code{ties} argument which is passed to \code{\link{estimate_p_value}}. For
#' options see \code{\link{estimate_p_value}}.
#'
#' The 100(1-alpha)\% global envelope is given by the 'k_alpha'th lower and
#' upper envelope. For details see Myllymäki et al. (2017).
#'
#' The above holds for p-value calculation if \code{lexo == FALSE} and then the test
#' corresponds to the rank envelope test by Myllymaki et. al (2013). If \code{lexo == TRUE},
#' then all the pointwise ranks are used to rank the curves by rank count ordering (Myllymäki et al., 2017)
#' and the single p-value in \code{p} is the p-value based on the rank count ordering.
#'
#' The rank count ordering test allows in principle a lower number of simulations to be used,
#' but then the test may no longer be usable as a graphical test.
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Seijo, H., Grabarnik, P. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. Global envelope tests for spatial point patterns. arXiv:1307.0239v4 [stat.ME]
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2013) Global envelope tests for spatial point patterns. arXiv:1307.0239v1 [stat.ME]
#'
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or an \code{\link[spatstat]{envelope}}
#'  object. If an envelope object is given, it must contain the summary
#'  functions from the simulated patterns which can be achieved by setting
#'  savefuns = TRUE when calling \code{\link[spatstat]{envelope}}.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param savedevs Logical. Should the global rank values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @param alternative A character string specifying the alternative hypothesis. Must be one of the following:
#'         "two.sided" (default), "less" or "greater".
#' @param lexo Logical, whether or not to use rank count ordering for calculation of the p-value. See details.
#' @param ties Ties method to be passed to \code{\link{estimate_p_value}}. Used to obtain
#' a point estimate for the p-value. The default point estimate is the mid-rank p-value.
#'
#' @return An object of class "envelope_test", "envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item data_curve = values of the test function for the data point pattern
#' \item lower = the lower envelope based on the simulated functions
#' \item upper = the upper envelope based on the simulated functions
#' \item central_curve = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Rank envelope test" for the rank envelope test)
#'   \item alternative = The alternative specified in the function call.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item p_interval = The p-value interval [p_liberal, p_conservative].
#'   \item ties = As the argument \code{ties}.
#'   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
#'   \item k = Global rank values (k[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type.
#' @export
#' @seealso \code{\link{random_labelling}}, \code{\link{plot.envelope_test}}
#' @examples
#'
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- unmark(spruces)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=2499, savefuns=TRUE, correction="translate")
#' # The rank envelope test
#' res <- rank_envelope(env)
#' # Plot the result.
#' # - The central curve is now obtained from env[['theo']], which is the
#' # value of the L-function under the null hypothesis (L(r) = r).
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=TRUE)
#'
#' ## Advanced use:
#' # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#' curve_set <- crop_curves(env, r_min = 1, r_max = 7)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # Do the rank envelope test
#' res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # 1) Perform simulations under the random labelling hypothesis and calculate
#' # the test function T(r) for the data pattern (mpp) and each simulation.
#' # The command below specifies that the test function is T(r) = \hat{L}_m(r),
#' # which is an estimator of the mark-weighted L function, L_m(r),
#' # with translational edge correction (default).
#' # The random_labelling function returns the centred functions \hat{L}_m(r)-T_0(r),
#' # where T_0(r) = \hat{L}(r) is the unmarked L function.
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' # 2) Do the rank envelope test
#' res <- rank_envelope(curve_set)
#' # 3) Plot the test result
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
#'
#' # Make the test using instead the test function T(r) = \hat{L}_mm(r);
#' # which is an estimator of the mark-weighted L function, L_mm(r),
#' # with translational edge correction (default).
#' curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- rank_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[mm](r)-L(r))))
#'
#' ## Goodness-of-fit test (typically conservative)
#' #-----------------------------------------------
#' pp <- unmark(spruces)
#' # Minimum distance between points in the pattern
#' min(nndist(pp))
#' # Fit a model
#' fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process
#'
#' \dontrun{
#' # Simulating Gibbs process by 'envelope' is slow, because it uses the MCMC algorithm
#' #env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE,
#'                  correction="none", r=seq(0, 4, length=500))
#'
#' # Using direct algorihm can be faster, because the perfect simulation is used here.
#' simulations <- NULL
#' for(j in 1:2499) {
#'    simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]),
#'                                  R = fittedmodel$interaction$par$hc,
#'                                  W = pp$window);
#'    if(j%%10==0) cat(j, "...", sep="")
#' }
#' env <- envelope(pp, simulate=simulations, fun="Jest", nsim=length(simulations),
#'                 savefuns=TRUE, correction="none", r=seq(0, 4, length=500))
#' curve_set <- crop_curves(env, r_min = 1, r_max = 3.5)
#' res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#' }
#'
rank_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE,
                          alternative=c("two.sided", "less", "greater"),
                          lexo=FALSE, ties) {
    # saving for attributes / plotting purposes
    lo.name <- "lower critical boundary for %s"
    hi.name <- "upper critical boundary for %s"
    switch(alternative,
            two.sided = {},
            less = {
                hi.name <- "infinite upper boundary"
            },
            greater = {
                lo.name <- "infinite lower boundary"
            })
    if(inherits(curve_set, 'envelope')) {
        fname <- attr(curve_set, "fname")
        labl <- attr(curve_set, "labl")
        desc <- attr(curve_set, "desc")
        desc[4] <- lo.name
        desc[5] <- hi.name
    }
    else {
        fname <- "T"
        labl <- c("r", "T[obs](r)", "T[0](r)", "T[lo](r)", "T[hi](r)")
        desc <- c("distance argument r",
                  "observed value of %s for data pattern",
                  "central curve under the null hypothesis",
                  lo.name, hi.name)
    }

    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")
    alternative <- match.arg(alternative)

    # The type of the p-value
    if(missing(ties))
        ties <- p_value_ties_default()
    else if(lexo) cat("The argument ties ignored, because lexo = TRUE. \n")
    if(lexo) ties <- "lexical"

    # data_curve = the vector of test function values for data
    # sim_curves = matrix where each row contains test function values of a simulation under null hypothesis
    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- length(curve_set$r)
    # Define the central curve T_0
    T_0 <- get_T_0(curve_set)

    data_and_sim_curves <- rbind(data_curve, sim_curves)
    loranks <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
    hiranks <- Nsim+2-loranks
    # k:
    switch(alternative,
           "two.sided" = {
               allranks <- pmin(loranks, hiranks)
           },
           "less" = {
               allranks <- loranks
           },
           "greater" = {
               allranks <- hiranks
           })

    distance <- apply(allranks, MARGIN=1, FUN=min)
    u <- -distance
    #-- p-interval
    p_low <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='liberal')
    p_upp <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='conservative')

    #-- p-value
    if(!lexo) {
        p <- estimate_p_value(x=u[1], sim_vec=u[-1], ties=ties)
    }
    else { # rank the curves by lexical ordering
        # order ranks within each curve
        sortranks <- apply(allranks, 1, sort) # curves now represented as columns
        lexo_values <- do.call("order", split(sortranks, row(sortranks)))

        # find ties
        sorted <- sortranks[ ,lexo_values]
        dupp <- duplicated(split(sorted, col(sorted)))
        tied <- dupp | c(dupp[-1], FALSE)

        # replace ranks of tied values by mean ranks
        # (maybe a little bit awkward, but isntitcool)
        tie.rle <- rle(tied)
        tie.end <- cumsum(tie.rle$lengths)
        tie.start <- cumsum(c(1,tie.rle$lengths))
        tie.start <- tie.start[-length(tie.start)]
        rank.rle <- tie.rle
        rank.rle$values <- (tie.start + tie.end)/2
        tieranks <- inverse.rle(rank.rle)
        newranks <- 1:(Nsim+1)
        newranks[tied] <- tieranks[tied]

        distance_lexo <- newranks[order(lexo_values)]
        #-- calculate the p-value
        u_lexo <- -distance_lexo
        p <- estimate_p_value(x=u_lexo[1], sim_vec=u_lexo[-1])
    }

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance, decreasing=TRUE)
    kalpha <- distancesorted[floor((1-alpha)*(Nsim+1))]
    LB <- array(0, nr);
    UB <- array(0, nr);

    for(i in 1:nr){
        Hod <- sort(data_and_sim_curves[,i])
        LB[i]<- Hod[kalpha];
        UB[i]<- Hod[Nsim+1-kalpha+1];
    }

    res <- structure(list(r=curve_set[['r']], obs=data_curve, theo=T_0, lo=LB, hi=UB),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Rank envelope test"
    attr(res, "alternative") <- alternative
    attr(res, "p") <- p
    attr(res, "p_interval") <- c(p_low, p_upp)
    attr(res, "ties") <- ties
    attr(res, "k_alpha") <- kalpha
    if(savedevs) attr(res, "k") <- distance
    # for fv
    attr(res, "fname") <- fname
    attr(res, "argu") <- "r"
    attr(res, "valu") <- "data_curve"
    attr(res, "ylab") <- "%s(r)"
    attr(res, "fmla") <- ". ~ r"
    attr(res, "alim") <- c(min(curve_set[['r']]), max(curve_set[['r']])) # FIXME, ?
    attr(res, "labl") <- labl
    attr(res, "desc") <- desc
    #attr(res, "unitname") <- "unit / units"
    attr(res, "shade") <- c("lower", "upper")
    attr(res, "call") <- match.call()
    res
}

#' Print method for the class 'envelope_test'
#' @usage \method{print}{envelope_test}(x, ...)
#'
#' @param x an 'envelope_test' object
#' @param ... Ignored.
#'
#' @method print envelope_test
#' @export
print.envelope_test <- function(x, ...) {
    cat(attr(x, "method"), "\n",
        " p-value of the test: ", attr(x, "p"), sep="")
    if(!is.null(attr(x, "ties"))) cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else cat("\n")
    if(!is.null(attr(x, "p_interval")))
        cat(" p-interval         : (", attr(x, "p_interval")[1], ", ", attr(x, "p_interval")[2],")\n", sep="")
}

#' Plot method for the class 'envelope_test'
#' @usage \method{plot}{envelope_test}(x, use_ggplot2=FALSE, base_size=15, dotplot=length(x$r)<10,
#'                                      main, ylim, xlab, ylab, ...)
#'
#' @param x an 'envelope_test' object
#' @param use_ggplot2 TRUE/FALSE, If TRUE, then a plot with a coloured envelope ribbon is provided. Requires R library ggplot2.
#' @param base_size Base font size, to be passed to theme style when \code{use_ggplot2 = TRUE}.
#' @param dotplot Logical. If TRUE, then instead of envelopes a dot plot is done.
#' Suitable for low dimensional test vectors. Only applicable if \code{use_ggplot2} is FALSE.
#' Default: TRUE if the dimension is less than 10, FALSE otherwise.
#' @param main See \code{\link{plot.default}}. A sensible default exists.
#' @param ylim See \code{\link{plot.default}}. A sensible default exists.
#' @param xlab See \code{\link{plot.default}}. A sensible default exists.
#' @param ylab See \code{\link{plot.default}}. A sensible default exists.
#' @param ... Additional parameters to be passed to \code{\link{env_basic_plot}}, \code{\link{dotplot}}
#' (if dotplot=TRUE) or \code{\link{env_ggplot}} (if use_ggplot2=TRUE).
#'
#' @method plot envelope_test
#' @export
#' @seealso \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}
plot.envelope_test <- function(x, use_ggplot2=FALSE, base_size=15, dotplot=length(x$r)<10,
        main, ylim, xlab, ylab, ...) {
    if(missing('main')) main <- env_main_default(x)
    if(missing('ylim')) ylim <- env_ylim_default(x, use_ggplot2)
    if(missing('xlab')) xlab <- expression(italic(r))
    if(missing('ylab')) ylab <- expression(italic(T(r)))

    if(use_ggplot2 & attr(x, "alternative") == "two.sided") {
        env_ggplot(x, base_size, main, ylim, xlab, ylab, ...)
    }
    else {
        if(use_ggplot2) cat("The use_ggplot2 option is valid only for the alternative \'two.sided\'. use_ggplot2 ignored.\n")
        if(dotplot) {
            warning("The plot style \'dotplot'\ does not search for combined tests.\n")
            env_dotplot(x, main, ylim, xlab, ylab, ...)
        }
        else {
            env_basic_plot(x, main, ylim, xlab, ylab, ...)
        }
    }
}

#' Studentised envelope test
#'
#' The studentised envelope test, which takes into account the unequal
#' variances of the test function T(r) for different distances r.
#'
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028 [stat.ME]
#'
#' Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P. (2013). Global envelope tests for spatial point patterns. arXiv:1307.0239 [stat.ME]
#'
#' @inheritParams rank_envelope
#' @param savedevs Logical. Should the deviation values u_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @param ... Additional parameters passed to \code{\link{estimate_p_value}} to obtain a point estimate
#' for the p-value. The default point estimate is the mid-rank p-value. The choice should not affect the
#' result, since no ties are expected to occur.
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r = Distances for which the test was made.
#'   \item method = The name of the envelope test.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item central_curve = If the curve_set (or envelope object) contains a component 'theo',
#'         then this function is used as the central curve and returned in this component.
#'         Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'   \item data_curve = The test function for the data.
#'   \item lower = The lower envelope.
#'   \item upper = The upper envelope.
#'   \item call = The call of the function.
#' }
#' @export
#' @importFrom stats sd
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The studentised envelope test
#' res <- st_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=TRUE)
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The studentised envelope test
#' res <- st_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- st_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
st_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, ...) {

    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    sdX <- as.vector(apply(curve_set[['sim_m']], MARGIN=1, FUN=stats::sd))

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    scaled_curve_set <- weigh_curves(curve_set, divisor_to_coeff(sdX))
    #devs <- deviation(scaled_curve_set, measure = 'max', scaling='qdir')
    # u_1
    distance[1] <- max(abs(scaled_curve_set$obs))
    # u_2, ..., u_{s+1}
    distance[2:(Nsim+1)] <- apply(abs(scaled_curve_set[['sim_m']]), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1], ...)

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance);
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))];
    LB <- T_0 - talpha*sdX;
    UB <- T_0 + talpha*sdX;

    res <- structure(list(r=curve_set[['r']], data_curve=data_curve, lower=LB, upper=UB, central_curve=T_0),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Studentised envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    attr(res, "call") <- match.call()
    res
}

#' Directional quantile envelope test
#'
#' The directional quantile envelope test, which takes into account the unequal 
#' variances of the test function T(r) for different distances r and is also 
#' protected against asymmetry of T(r).
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028 [stat.ME]
#'
#' Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P. (2013). Global envelope tests for spatial point patterns. arXiv:1307.0239 [stat.ME]
#'
#' @inheritParams st_envelope
#' @param probs A two-element vector containing the lower and upper
#'   quantiles for the envelope, in that order and on the interval [0, 1].
#'   The default values are 0.025 and 0.975.
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r = Distances for which the test was made.
#'   \item method = The name of the envelope test.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item central_curve = If the curve_set (or envelope object) contains a component 'theo',
#'         then this function is used as the central curve and returned in this component.
#'         Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'   \item data_curve = The test function for the data.
#'   \item lower = The lower envelope.
#'   \item upper = The upper envelope.
#'   \item call = The call of the function.
#' }
#' @export
#' @importFrom stats quantile
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The directional quantile envelope test
#' res <- qdir_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=TRUE)
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The directional quantile envelope test
#' res <- qdir_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- qdir_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
qdir_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, probs = c(0.025, 0.975), ...) {

    curve_set <- convert_envelope(curve_set)
    check_probs(probs)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    # calculate quantiles for residual curve_set (i.e. for sim_curves - T_0)
    quant_m <- apply(curve_set[['sim_m']], 1, stats::quantile, probs = probs)
    abs_coeff <- divisor_to_coeff(abs(quant_m))
    lower_coeff <- abs_coeff[1, , drop = TRUE]
    upper_coeff <- abs_coeff[2, , drop = TRUE]

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    # u_1
    scaled_residuals <- weigh_both_sides(curve_set[['obs']], upper_coeff, lower_coeff)
    distance[1] <- max(abs(scaled_residuals))
    # u_2, ..., u_{s+1}
    sim_scaled_residuals <- weigh_both_sides(curve_set[['sim_m']], upper_coeff, lower_coeff)
    distance[2:(Nsim+1)] <- apply(abs(sim_scaled_residuals), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1], ...)

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance)
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))]
    LB <- T_0 - talpha*abs(quant_m[1,])
    UB <- T_0 + talpha*abs(quant_m[2,])

    res <- structure(list(r=curve_set[['r']], data_curve=data_curve, lower=LB, upper=UB, central_curve=T_0),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Directional quantile envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    attr(res, "call") <- match.call()
    res
}

#' Unscaled envelope test
#'
#' The unscaled envelope test, which leads to envelopes with constant width
#' over the distances r. It corresponds to the classical maximum deviation test
#' without scaling.
#'
#'
#' This test suffers from unequal variance of T(r) over the distances r and from
#' the asymmetry of distribution of T(r). We recommend to use the rank_envelope
#' (if number of simulations close to 5000 can be afforded) or st_envelope/qdir_envelope
#' (if large number of simulations cannot be afforded) instead.
#'
#' @references
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @inheritParams st_envelope
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r = Distances for which the test was made.
#'   \item method = The name of the envelope test.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern). Returned only if savedevs = TRUE.
#'   \item central_curve = If the curve_set (or envelope object) contains a component 'theo',
#'         then this function is used as the central curve and returned in this component.
#'         Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'   \item data_curve = The test function for the data.
#'   \item lower = The lower envelope.
#'   \item upper = The upper envelope.
#'   \item call = The call of the function.
#' }
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The studentised envelope test
#' res <- unscaled_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=TRUE)
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 8)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The studentised envelope test
#' res <- unscaled_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#'
#' ## Random labeling test
#' #----------------------
#' # requires library 'marksummary'
#' mpp <- spruces
#' # Use the test function T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
#' res <- unscaled_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
unscaled_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, ...) {

    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- dim(sim_curves)[2]

    # Define T_0 and residual curve_set
    T_0 <- get_T_0(curve_set)
    curve_set <- residual(curve_set, use_theo = TRUE)

    # Calculate deviation measures
    distance <- array(0, Nsim+1);
    # u_1
    distance[1] <- max(abs(curve_set$obs))
    # u_2, ..., u_{s+1}
    distance[2:(Nsim+1)] <- apply(abs(curve_set[['sim_m']]), 2, max)

    #-- calculate the p-value
    p <- estimate_p_value(x=distance[1], sim_vec=distance[-1], ...)

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance);
    talpha <- distancesorted[floor((1-alpha)*(Nsim+1))];
    LB <- T_0 - talpha;
    UB <- T_0 + talpha;

    res <- structure(list(r=curve_set[['r']], data_curve=data_curve, lower=LB, upper=UB, central_curve=T_0),
                     class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Unscaled envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    if(savedevs) attr(res, "u") <- distance
    attr(res, "call") <- match.call()
    res
}


#' Approximative normal envelope test
#'
#' Approximative normal envelope test
#'
#'
#' The normal envelope test is a parametric envelope test.
#' If simulation from the null model is extremely tedious, for some test functions T(r)
#' it is possible to use a normal approximation as suggested by Mrkvicka (2009) in a study
#' of random closed sets.
#' The basis of the test is the approximation of the test function T(r), r in I=[r_min, r_max],
#' by a random vector (T(r_1), ...,T(r_m))', where m is the number of distances, that follows
#' multivariate normal distribution with mean mu and variance matrix Sigma which are estimated
#' from T_j(r), j=2, ...,s+1. This test employes T_j(r), j=2, ...,s+1, only for estimating mu
#' and Sigma, and for this s=100 simulations is enough to reach needed accuracy.
#'
#' This normal envelope test is not a Monte Carlo test. It employes the same envelopes as the
#' studentised envelope test (see \code{\link{st_envelope}}), i.e. the semiparametric kth lower
#' and upper envelopes
#'
#' T^u_{low}(r)= T_0(r) - u sqrt(var_0(T(r)))
#' and
#' T^u_{upp}(r)= T_0(r) + u sqrt(var_0(T(r))),
#'
#' but the p-value and the 100(1-alpha) percent global envelope are calculated based
#' on simulations from a multivariate normal distribution.
#'
#'
#' @references Mrkvička, T. (2009). On testing of general random closed set model hypothesis. Kybernetika 45, 293-308.
#' @inheritParams st_envelope
#' @param n_norm Number of simulations drawn from the multivariate normal distribution (dimension = number of distances r).
#' @export
#' @importFrom stats var
#' @importFrom stats sd
#' @importFrom mvtnorm rmvnorm
#' @examples
#' require(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=99, savefuns=TRUE,
#'                 correction="translate", r=seq(0,8,length=50))
#' curve_set <- residual(env, use_theo = TRUE)
#' system.time( res <- normal_envelope(curve_set, n_norm=200000) )
#' plot(res)
normal_envelope <- function(curve_set, alpha=0.05, n_norm=200000, ...) {
    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    n <-length(data_curve);
    EX <- colMeans(sim_curves, na.rm = TRUE);
    varX <- stats::var(sim_curves, na.rm = TRUE);

    #-- simulate from the normal distribution
    simnorm <- mvtnorm::rmvnorm(n=n_norm, mean = EX, sigma = varX, method=c('svd'));

    sdX <- as.vector(apply(sim_curves, MARGIN=2, FUN=stats::sd))
    distance <- array(0, n_norm);
    for(j in 1:n_norm) {
        ttt <- abs(simnorm[j,]-EX)/sdX
        ttt[!is.finite(ttt)] <- 0
        distance[j] <- max(ttt);
    }

    distancesorted <- sort(distance);

    ttt<-abs(data_curve-EX)/sdX;
    ttt[!is.finite(ttt)] <- 0
    tmaxd <- max(ttt)

    #-- calculate the p-value
    p <- estimate_p_value(x=tmaxd, sim_vec=distance, ...)
    #    p <- 1;
    #    for(i in 1:n_norm) {
    #        if(distance[i]>tmaxd) p<-p+1;
    #    }
    #    p <- p/(n_norm+1);

    #-- calculate the 100(1-alpha)% global envelope
    talpha <- distancesorted[floor((1-alpha)*n_norm)];
    LB <- EX - talpha*sdX
    UB <- EX + talpha*sdX

    res <- structure(list(r=curve_set[['r']], data_curve=data_curve, lower=LB, upper=UB, central_curve=EX),
            class = c("envelope_test", "envelope", "fv", "data.frame"))
    attr(res, "method") <- "Approximative normal envelope test"
    attr(res, "alternative") <- "two.sided"
    attr(res, "p") <- p
    attr(res, "u_alpha") <- talpha
    attr(res, "call") <- match.call()
    res
}
