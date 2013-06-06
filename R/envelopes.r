#' The rank envelope test
#'
#' @inheritParams convert_envelope
#' @param alpha
#' @param ... Additional parameters passed to \code{\link{estimate_p_value}}.
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' # Test of complete spatial randomness (CSR)
#' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=4999, savefuns=TRUE, correction="translate")
#' # Create a curve set on the interval of distances [r_min, r_max]
#' curve_set <- crop_curves(env, r_min = 1, r_max = 7)
#' # The rank envelope test
#' res <- rank_envelope(curve_set)
#' plot(res)
#'
#' # For visualisation, take the L(r)-r function
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' # The result of the rank envelope test stays the same
#' res <- rank_envelope(curve_set); plot(res)
rank_envelope <- function(curve_set, alpha=0.05, ...) {
    # data_curve = the vector of L-function values for data
    # sim_curves = matrix where each row contains L function values of a simulation under null hypothesis
    # Nsim  = number of simulations, that is nrow(sim_curves) --> This argument not needed
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
    p_low <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='liberal')
    p_upp <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='conservative')

    #-- calculate the 100(1-alpha)% envelopes
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
                central_curve=MX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
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
                    p <- p + geom_hline(yintercept=0, color='grey30', linetype='dashed', size=0.1)
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

#' Variance stabilized envelope test
#'
#' @inheritParams rank_envelope
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
#' curve_set <- crop_curves(env, r_min = 0, r_max = 8)
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' res <- sd_envelope(curve_set)
#' plot(res)
#' plot(res, xlab=expression(italic(r)), ylab=expression(italic(L(r)-r)))
sd_envelope <- function(curve_set, alpha=0.05, ...) {

    curve_set <- convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    n <- dim(sim_curves)[2]
    EX <- colMeans(sim_curves);
    sdX <- as.vector(apply(sim_curves, MARGIN=2, FUN=sd))

    #    for(i in 1:n) {
    #        if(sdX[i]==0) {
    #            sdX[i]<-0.00000001
    #        }
    #    }

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

    #-- calculate the 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*(Nsim+1))];
    LB <- EX - talpha*sdX;
    UB <- EX + talpha*sdX;

    res <- list(r=curve_set[['r']], method="Variance stabilized envelope test", p=p,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    class(res) <- "envelope_test"
    res
}

#' Asymmetric quantile envelope test
#'
#' @inheritParams rank_envelope
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
#' curve_set <- crop_curves(env, r_min = 0, r_max = 8)
#' curve_set <- residual(curve_set, use_theo = TRUE)
#' res <- as_quantile_envelope(curve_set)
#' plot(res)
as_quantile_envelope <- function(curve_set, alpha=0.05, ...) {

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

    #-- calculate the 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*(Nsim+1))];
    LB <- EX - talpha*abs(QQ[1,]-EX);
    UB <- EX + talpha*abs(QQ[2,]-EX);

    res <- list(r=curve_set[['r']], method="Asymmetric quantile envelope test", p=p,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    class(res) <- "envelope_test"
    res
}



#' Approximative normal envelope test
#'
#' @inheritParams rank_envelope
#' @param n_norm Number of simulations drawn from the multivariate normal distribution (dimension = number of distances r).
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=99, savefuns=TRUE, r=seq(0,8,length=50))
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
    #    sdX <- array(0,n);
    #    for(i in 1:n){sdX[i] <- sqrt(varX[i,i])}
    #    for(i in 1:n){if(sdX[i]==0){sdX[i]<-0.00000001}}
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

    #-- calculate the 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*n_norm)];
    LB <- EX - talpha*sdX
    UB <- EX + talpha*sdX

    res <- list(r=curve_set[['r']], method="Approximative normal envelope test", p=p,
                central_curve=EX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    class(res) <- "envelope_test"
    res
}
