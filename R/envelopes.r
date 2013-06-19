#' The rank envelope test
#'
#' The rank envelope test, p-value and exact simultaneous 100(1-alpha) percent envelopes
#'
#'
#' The rank envelope test is a completely non-parametric test, which provides a p-value
#' interval given by the most liberal and the most conservative p-value estimate and
#' the simultaneous 100(1-alpha) percent envelopes for the chosen test function T(r) on
#' the chosen interval of distances.
#'
#' @references Myllymäki, M., Mrkvička, T., Seijo, H., Grabarnik, P. (2013). Global envelope tests for spatial point patterns.
#'
#' @inheritParams convert_envelope
#' @param alpha The significance level. Simultaneous 100(1-alpha) percent envelopes will be calculated.
#' @param savedevs Logical. Should the global rank values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @param ... Additional parameters passed to \code{\link{estimate_p_value}} to obtain a point estimate for the p-value.
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r Distances for which the test was made.
#'   \item method The envelope method.
#'   \item p A point estimate for the p-value.
#'   \item p_interval The p-value interval [p_liberal, p_conservative].
#'   \item central_curve The mean test function (median) calculated from simulations.
#'   \item data_curve The test function for the data.
#'   \item lower The lower envelope.
#'   \item upper The upper envelope.
#'   \item call The call of the function.
#' }
#' @export
#' @seealso \code{\link{random_labeling}}
#' @examples
#'
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- unmark(spruces)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=4999, savefuns=TRUE, correction="translate")
#' # The rank envelope test
#' res <- rank_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=TRUE)
#'
#' ## Advanced use:
#' # Create a curve set, choosing the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 7)
#' # For better visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The rank envelope test
#' res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
#'
#' ## Random labeling test
#' #----------------------
#' mpp <- spruces
#' # T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
#' res <- rank_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
#' # T(r) = \hat{L}_mm(r), an estimator of the L_mm(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=4999, r_min=0, r_max=9.5)
#' res <- rank_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[mm](r)-L(r))))
#'
#' ## Goodness-of-fit test
#' #----------------------
#' # FIXME What is a reasonable model to be fitted to this data?
#' pp <- unmark(spruces)
#' # Minimum distance between points in the pattern
#' min(nndist(pp))
#' # Fit a model
#' ## fittedmodel <- ppm(pp, interaction=Strauss(r=4)) # Strauss process
#' fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process
#'
#' \dontrun{
#' ## Simulating Strauss process by 'envelope' is slow
#' #env <- envelope(fittedmodel, fun="Lest", nsim=999, savefuns=TRUE, correction="none", r=seq(0, 9.5, length=500))
#'
#' simulations <- NULL
#' for(j in 1:999) {
#'    ##simulations[[j]] <- rStrauss(beta=exp(fittedmodel$coef[1]), gamma=exp(fittedmodel$coef[2]), R=fittedmodel$interaction$par$r, W=pp$window);
#'    simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]), R = fittedmodel$interaction$par$hc, W = pp$window);
#'    if(j%%10==0) cat(j, "...", sep="")
#' }
#' env <- envelope(pp, simulate=simulations, fun="Lest", nsim=length(simulations), savefuns=TRUE, correction="none", r=seq(0, 9.5, length=500))
#' res <- rank_envelope(env); plot(res, use_ggplot2=TRUE)
#' }
#'
rank_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, ...) {
    # data_curve = the vector of L-function values for data
    # sim_curves = matrix where each row contains L function values of a simulation under null hypothesis
    # alpha = the chosen significance level of the test

    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    n <- length(curve_set$r)
    MX <- apply(sim_curves, MARGIN=2, FUN=median)
    data_and_sim_curves <- rbind(data_curve, sim_curves)
    RR <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
    Rmin <- apply(RR, MARGIN=1, FUN=min)
    Rmax <- Nsim+2-apply(RR, MARGIN=1, FUN=max)
    RmRm <- rbind(Rmin,Rmax)
    # k:
    distance <- apply(RmRm, MARGIN=2, FUN=min)

    #-- calculate the p-value
    u <- -distance
    p <- estimate_p_value(obs=u[1], sim_vec=u[-1], ...)
    # p-interval
    p_low <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='liberal')
    p_upp <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='conservative')

    #-- calculate the simultaneous 100(1-alpha)% envelopes
    distancesorted <- sort(distance, decreasing=TRUE)
    kalpha <- distancesorted[round((1-alpha)*(Nsim+1))]
    LB <- array(0, n);
    UB <- array(0, n);

    for(i in 1:n){
        Hod <- sort(data_and_sim_curves[,i])
        LB[i]<- Hod[kalpha];
        UB[i]<- Hod[Nsim+1-kalpha+1];
    }

    res <- list(r=curve_set[['r']], method="Rank envelope test",
                p=p, p_interval=c(p_low,p_upp),
                k_alpha=kalpha,
                central_curve=MX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    if(savedevs) res$k <- distance
    class(res) <- "envelope_test"
    res
}

#' Print method for the class 'envelope_test'
#' @usage \method{print}{envelope_test}(x)
#'
#' @param x an 'envelope_test' object
#'
#' @method print envelope_test
#' @export
print.envelope_test <- function(x, ...) {
    with(x, cat(method, "\n",
                "p-value of the test:", p, "\n"))
    with(x, if(exists('p_interval'))
            cat(" p-interval         : (", p_interval[1], ", ", p_interval[2],")\n", sep=""))
}

#' Plot method for the class 'envelope_test'
#' @usage \method{plot}{envelope_test}(x)
#'
#' @param x an 'envelope_test' object
#' @param use_ggplot2 TRUE/FALSE, If TRUE, then a plot with a coloured envelope ribbon is provided. Requires R library ggplot2.
#' @param main See \code{\link{plot.default}}. A sensible default exists.
#' @param ylim See \code{\link{plot.default}}. A sensible default exists.
#' @param xlab See \code{\link{plot.default}}. A sensible default exists.
#' @param ylab See \code{\link{plot.default}}. A sensible default exists.
#' @param ... Additional parameters to be passed to plot (if use_ggplot2=FALSE).
#'
#' @method plot envelope_test
#' @export
#' @seealso \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}
plot.envelope_test <- function(x, use_ggplot2=FALSE, main, ylim, xlab, ylab, ...) {
    if(missing('main')) {
        if(with(x, exists('p_interval')))
            main <- paste(x$method, ": p-interval = (",
                    round(x$p_interval[1],3),", ", round(x$p_interval[2],3), ")", sep="")
        else
            main <- paste(x$method, ": p = ", round(x$p,3), sep="")
    }
    if(missing('ylim')) {
        ylim <- c(min(x$data_curve,x$lower,x$upper,x$central_curve),
                  max(x$data_curve,x$lower,x$upper,x$central_curve))
    }
    if(missing('xlab')) xlab <- expression(italic(r))
    if(missing('ylab')) ylab <- expression(italic(T(r)))

    if(use_ggplot2) {
        require(ggplot2)
        linetype.values <- c('solid', 'dashed')
        with(x, {
                    df <- data.frame(r = rep(r, times=2),
                                     curves = c(data_curve, central_curve),
                                     type = factor(rep(c("Data function", "Central function"), each=length(r)), levels=c("Data function", "Central function")),
                                     lower = rep(lower, times=2),
                                     upper = rep(upper, times=2),
                                     main = factor(rep(main, times=length(r)))
                                     )
                    p <- (ggplot(df, aes_string(x='r', y='curves', group='type', linetype='type'))
                                + geom_line(aes_string(linetype='type'))
                                + geom_ribbon(aes(ymin=lower, ymax=upper), fill='grey59',
                                        alpha=1)
                                + geom_line(aes(y=curves))
                                + facet_grid('~ main', scales='free')
                                + scale_x_continuous(name=xlab)
                                + scale_y_continuous(name=ylab)
                                + scale_linetype_manual(values=linetype.values, name='')
                                + ThemePlain()
                                )
                    #p <- p + geom_hline(yintercept=0, color='grey30', linetype='dashed', size=0.1)
                    print(p)
                }
            )
    }
    else {
        with(x, {
                    plot(r, data_curve, ylim=ylim, main=main, xlab=xlab, ylab=ylab,
                            type="l", lty=1, lwd=2, ...)
                    lines(r, lower, lty=2)
                    lines(r, upper, lty=2)
                    lines(r, central_curve, lty=1)
                }
        )
    }
}

#' Studentised envelope test
#'
#' Studentised envelope test
#'
#'
#' The studentised envelope test is a semiparametric envelope test,
#' which takes into account the unequal variances of the test function
#' T(r) for different distances r.
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
#'
#' Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P. (2013). Global envelope tests for spatial point patterns.
#'
#' @inheritParams rank_envelope
#' @param savedevs Logical. Should the global deviation values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r Distances for which the test was made.
#'   \item method The envelope method.
#'   \item p A point estimate for the p-value.
#'   \item central_curve The mean test function (median) calculated from simulations.
#'   \item data_curve The test function for the data.
#'   \item lower The lower envelope.
#'   \item upper The upper envelope.
#'   \item call The call of the function.
#' }
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' library(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The studentised envelope test
#' res <- st_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=T)
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
#' mpp <- spruces
#' # T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
#' res <- st_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
st_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, ...) {

    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    n <- dim(sim_curves)[2]
    EX <- colMeans(sim_curves);
    sdX <- as.vector(apply(sim_curves, MARGIN=2, FUN=sd))

    distance <- array(0, Nsim+1);
    # data
    ttt <- abs(data_curve-EX)/sdX;
    ttt[!is.finite(ttt)] <- 0
    distance[1] <- max(ttt)
    # simulations
    for(j in 1:Nsim) {
        ttt <- abs(sim_curves[j,]-EX)/sdX
        ttt[!is.finite(ttt)] <- 0
        distance[j+1] <- max(ttt);
    }

    distancesorted <- sort(distance);

    #-- calculate the p-value
    p <- estimate_p_value(obs=distance[1], sim_vec=distance[-1], ...)

    #-- calculate the simultaneous 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*(Nsim+1))];
    LB <- EX - talpha*sdX;
    UB <- EX + talpha*sdX;

    res <- list(r=curve_set[['r']], method="Studentised envelope test", p=p,
                u_alpha=talpha,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    if(savedevs) res$u <- distance
    class(res) <- "envelope_test"
    res
}

#' Directional quantile envelope test
#'
#' Directional quantile envelope test
#'
#'
#' The directional quantile envelope test is a semiparametric envelope test,
#' which takes into account the unequal variances of the test function T(r)
#' for different distances r and is also protected against asymmetry of T(r).
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
#'
#' Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P. (2013). Global envelope tests for spatial point patterns.
#'
#' @inheritParams rank_envelope
#' @return An "envelope_test" object containing the following fields:
#' \itemize{
#'   \item r Distances for which the test was made.
#'   \item method The envelope method.
#'   \item p A point estimate for the p-value.
#'   \item central_curve The mean test function (median) calculated from simulations.
#'   \item data_curve The test function for the data.
#'   \item lower The lower envelope.
#'   \item upper The upper envelope.
#'   \item call The call of the function.
#' }
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' library(spatstat)
#' pp <- spruces
#' ## Test for complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The directional quantile envelope test
#' res <- qdir_envelope(env)
#' plot(res)
#' # or (requires R library ggplot2)
#' plot(res, use_ggplot2=T)
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
#' mpp <- spruces
#' # T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
#' res <- qdir_envelope(curve_set)
#' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
qdir_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE, ...) {

    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    n <- dim(sim_curves)[2]
    EX <- colMeans(sim_curves);
    QQ <- apply(sim_curves, MARGIN=2, FUN=quantile, probs = c(0.025,0.975))

    distance <- array(0, Nsim+1);
    tmaxd<-0;
    for(i in 1:length(EX)){
        if(data_curve[i]-EX[i]>0) {
            ttt <- (data_curve[i]-EX[i])/(QQ[2,i]-EX[i])
        }
        else {
            ttt <- (data_curve[i]-EX[i])/(QQ[1,i]-EX[i])
        }
        if(!is.finite(ttt)) ttt <- 0
        if(tmaxd<ttt) {
            tmaxd <- ttt
        }
    }
    distance[1] <- tmaxd
    for(j in 1:Nsim){
        tmax<-0;
        for(i in 1:length(EX)){
            if(sim_curves[j,i]-EX[i]>0) {
                ttt <- (sim_curves[j,i]-EX[i])/(QQ[2,i]-EX[i])
            }
            else {
                ttt <- (sim_curves[j,i]-EX[i])/(QQ[1,i]-EX[i])
            }
            if(!is.finite(ttt)) ttt <- 0
            if(tmax<ttt) {
                tmax <- ttt
            }
        }
        distance[j+1] <- tmax;
    }

    distancesorted <- sort(distance);

    #-- calculate the p-value
    p <- estimate_p_value(obs=distance[1], sim_vec=distance[-1], ...)

    #-- calculate the simultaneous 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*(Nsim+1))];
    LB <- EX - talpha*abs(QQ[1,]-EX);
    UB <- EX + talpha*abs(QQ[2,]-EX);

    res <- list(r=curve_set[['r']], method="Directional quantile envelope test", p=p,
                u_alpha=talpha,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    if(savedevs) res$u <- distance
    class(res) <- "envelope_test"
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
#' T^u_{low}(r)= T_0(r) - u sqrt{var_0(T(r))}
#' and
#' T^u_{upp}(r)= T_0(r) + u sqrt{var_0(T(r))},
#'
#' but the p-value and the simultaneous 100(1-alpha) percent envelopes are calculated based
#' on simulations from a multivariate normal distribution.
#'
#'
#' @references Mrkvička, T. (2009). On testing of general random closed set model hypothesis. Kybernetika 45, 293-308.
#' @inheritParams rank_envelope
#' @param n_norm Number of simulations drawn from the multivariate normal distribution (dimension = number of distances r).
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=99, savefuns=TRUE, correction="translate", r=seq(0,8,length=50))
#' curve_set <- residual(env, use_theo = TRUE)
#' system.time( res <- normal_envelope(curve_set, n_norm=200000) )
#' plot(res)
normal_envelope <- function(curve_set, alpha=0.05, n_norm=200000, ...) {
    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    n <-length(data_curve);
    EX <- colMeans(sim_curves, na.rm = TRUE);
    varX <- var(sim_curves, na.rm = TRUE);

    #-- simulate from the normal distribution
    simnorm <- rmvnorm(n=n_norm, mean = EX, sigma = varX, method=c('svd'));

    sdX <- as.vector(apply(sim_curves, MARGIN=2, FUN=sd))
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
    p <- estimate_p_value(obs=tmaxd, sim_vec=distance, ...)
    #    p <- 1;
    #    for(i in 1:n_norm) {
    #        if(distance[i]>tmaxd) p<-p+1;
    #    }
    #    p <- p/(n_norm+1);

    #-- calculate the simultaneous 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*n_norm)];
    LB <- EX - talpha*sdX
    UB <- EX + talpha*sdX

    res <- list(r=curve_set[['r']], method="Approximative normal envelope test", p=p,
                u_alpha = talpha,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    class(res) <- "envelope_test"
    res
}
