#' The rank envelope test
#'
#' @inheritParams convert_envelope
#' @param alpha
#' @param ... Additional parameters passed to \code{\link{estimate_p_value}}.
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
#' curve_set <- spptest:::crop_curves(env, r_min = 1, r_max = 7)
#' curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
#' res <- rank_envelope(curve_set)
#' plot(res)
rank_envelope <- function(curve_set, alpha=0.05, ...) {
    # Lm    = the vector of L-function values for data
    # stats = matrix where each row contains L function values of a simulation under null hypothesis
    # Nsim  = number of simulations, that is nrow(stats) --> This argument not needed
    # alpha = the chosen significance level of the test

    curve_set <- spptest:::convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    stats <- t(curve_set[['sim_m']])

    Nsim <- dim(stats)[1];
    n <- dim(stats)[2]
    MX <- apply(stats, MARGIN=2, FUN=median)
    statsLm <- rbind(stats,Lm);
    RR <- apply(statsLm, MARGIN=2, FUN=rank)
    Rmin <- apply(RR, MARGIN=1, FUN=min)
    Rmax <- apply(RR, MARGIN=1, FUN=max)
    Rmax <- Nsim+2-Rmax
    RmRm <- rbind(Rmin,Rmax)
    distance <- apply(RmRm, MARGIN=2, FUN=min)

    #-- calculate the p-value
    u <- -distance
    p <- spptest:::estimate_p_value(obs=u[Nsim+1], sim_vec=u[1:Nsim], ...)
    #    pm <- 1;pl <- 1;pu <- 1;
    #    for(i in 1:Nsim) {
    #        if (distance[i]<distance[Nsim+1]) {pm<-pm+1;pl<-pl+1;pu<-pu+1;}
    #        if (distance[i]==distance[Nsim+1]) {pm<-pm+1/2;pu<-pu+1;}
    #    }
    #    pm <- pm/(Nsim+1);
    #    pl <- pl/(Nsim+1);
    #    pu <- pu/(Nsim+1);

    #-- calculate the 100(1-alpha)% envelopes
    distancesorted <- sort(distance[-(Nsim+1)]);
    kalpha <- distancesorted[round((alpha)*Nsim)];
    LB <- array(0, n);
    UB <- array(0, n);

    for(i in 1:n){
        Hod <- sort(stats[,i])
        LB[i]<- Hod[kalpha];
        UB[i]<- Hod[Nsim-kalpha+1];
    }

    res <- list(r=curve_set[['r']], method="Rank envelope test", p=p,
                central_curve=MX, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    class(res) <- "envelope_test"
    res
}

print.envelope_test <- function(x, ...) {
    with(x, cat(method, "\n",
                "p-value of the test:", p, "\n"))
}

plot.envelope_test <- function(x, ...) {
    with(x, {
         plot(r, data_curve, ylim=c(min(data_curve,lower,upper,central_curve), max(data_curve,lower,upper,central_curve)), type="l", lty=1, lwd=2,
              main=method, ...) #, xlab=expression(r), ylab=expression(S(r))
         lines(r, lower, lty=2)
         lines(r, upper, lty=2)
         lines(r, central_curve, lty=1)
                })
}

#' Variance stabilized envelope test
#'
#' @inheritParams rank_envelope
#' @export
#' @examples
#' library(spatstat)
#' pp <- spruces
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
#' curve_set <- spptest:::crop_curves(env, r_min = 0, r_max = 8)
#' curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
#' res <- sd_envelope(curve_set)
#' plot(res)
#' plot(res, xlab=expression(italic(r)), ylab=expression(italic(L(r)-r)))
sd_envelope <- function(curve_set, alpha=0.05, ...) {

    curve_set <- spptest:::convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    stats <- t(curve_set[['sim_m']])

    Nsim <- dim(stats)[1];
    n <- dim(stats)[2]
    EX <- colMeans(stats);
    sdX <- as.vector(apply(stats, MARGIN=2, FUN=sd))

    #    for(i in 1:n) {
    #        if(sdX[i]==0) {
    #            sdX[i]<-0.00000001
    #        }
    #    }

    distance <- array(0, Nsim);
    for(j in 1:Nsim) {
        ttt <- abs(stats[j,]-EX)/sdX
        ttt[!is.finite(ttt)] <- 0
        distance[j] <- max(ttt);
    }

    distancesorted <- sort(distance);

    ttt <- abs(data_curve-EX)/sdX;
    ttt[!is.finite(ttt)] <- 0
    tmaxd <- max(ttt)

    #-- calculate the p-value
    p <- spptest:::estimate_p_value(obs=tmaxd, sim_vec=distance, ...)

    #-- calculate the 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*Nsim)];
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
#' curve_set <- spptest:::crop_curves(env, r_min = 0, r_max = 8)
#' curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
#' res <- as_quantile_envelope(curve_set)
#' plot(res)
as_quantile_envelope <- function(curve_set, alpha=0.05, ...) {

    curve_set <- spptest:::convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    stats <- t(curve_set[['sim_m']])

    Nsim <- dim(stats)[1];
    n <- dim(stats)[2]
    EX <- colMeans(stats);
    QQ <- apply(stats,MARGIN=2,FUN=quantile, probs = c(0.025,0.975))

    distance <- array(0, Nsim);
    for(j in 1:Nsim){
        tmax<-0;
        for(i in 1:length(EX)){
            if(stats[j,i]-EX[i]>0) {
                ttt <- (stats[j,i]-EX[i])/(QQ[2,i]-EX[i])
            }
            else {
                ttt <- (stats[j,i]-EX[i])/(QQ[1,i]-EX[i])
            }
            if(!is.finite(ttt)) ttt <- 0
            if(tmax<ttt) {
                tmax<-ttt
            }
        }
        distance[j] <- tmax;
    }

    distancesorted <- sort(distance);

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
            tmaxd<-ttt
        }
    }

    #-- calculate the p-value
    p <- spptest:::estimate_p_value(obs=tmaxd, sim_vec=distance, ...)

    #-- calculate the 100(1-alpha)% envelopes
    talpha <- distancesorted[round((1-alpha)*Nsim)];
    LB <- EX - talpha*(EX-QQ[1,]);
    UB <- EX + talpha*(QQ[2,]-EX);

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
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
#' curve_set <- spptest:::crop_curves(env, r_min = 0, r_max = 8)
#' curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
#' res <- normal_envelope(curve_set, n_norm=20000)
#' plot(res)
normal_envelope <- function(curve_set, alpha=0.05, n_norm=200000, ...) {
    curve_set <- spptest:::convert_envelope(curve_set)

    data_curve <- curve_set[['obs']]
    stats <- t(curve_set[['sim_m']])

    n <-length(data_curve);
    EX <- colMeans(stats, na.rm = TRUE);
    varX <- var(stats, na.rm = TRUE);

    #-- simulate from the normal distribution
    simnorm <- rmvnorm(n=n_norm, mean = EX, sigma = varX, method=c('svd'));

    sdX <- as.vector(apply(stats, MARGIN=2, FUN=sd))
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
    p <- spptest:::estimate_p_value(obs=tmaxd, sim_vec=distance, ...)
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
