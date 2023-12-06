# Matching function of Vilodomat et al. (2014) provided at https://doi.org/10.1111/biom.12139
# Tiny changes only.
#' @importFrom stats rnorm
#' @importFrom stats fitted
matching <- function( X.randomized, long, lat, Delta, target_variog, prctile, ids, maxk ) {
  variog.X.delta <- vector(mode = "list", length = length(Delta))
  linear.fit <- variog.X.delta
  hat.X.delta <- variog.X.delta
  resid.sum.squares <- rep(0, length(Delta))
  for(k in seq_along(Delta)) {
    # smooth X.randomized using locfit:
    fit <- locfit::locfit(X.randomized ~ locfit::lp(long, lat, nn = Delta[k], deg = 0), kern = "gauss", maxk = maxk)
    X.delta <- fitted(fit)
    # variogram of X.delta:
    variog.X.delta[[k]] <- geoR::variog(data = X.delta[ids], coords = cbind(long[ids],lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)

    # linear regression between the target and X.delta variograms:
    linear.fit[[k]] <- lm(target_variog$v ~ 1 + variog.X.delta[[k]]$v)
    # least square estimates:
    bet.hat <- as.numeric(linear.fit[[k]]$coefficients)
    # transformed X.delta:
    hat.X.delta[[k]] <- X.delta * sqrt(abs(bet.hat[2])) + rnorm(length(X.delta)) * sqrt(abs(bet.hat[1]))
    variog.hat.X.delta <- geoR::variog(data = hat.X.delta[[k]][ids], coords = cbind(long[ids],lat[ids]), option = "bin", max.dist = prctile, messages = FALSE)

    # sum of squares of the residuals:
    resid.sum.squares[k] <- sum((variog.hat.X.delta$v-target_variog$v) ^ 2)
  }
  # delta that minimizes the residual sum of squares:
  delta.star.id <- which.min(resid.sum.squares)
  hat.X.delta.star <- hat.X.delta[[delta.star.id]]
  return(list(residus = resid.sum.squares, delta.star.id = delta.star.id, hat.X.delta.star = hat.X.delta.star))
}

#' The test of local correlations
#'
#' The test of local correlations using Vilodomat et al. (2014) procedure for resamples
#' and the FDR envelope of Mrkvička and Myllymäki (2023).
#'
#'
#' The code is a modification of the supporting information code of Vilodomat et al. (2014)
#' available at https://doi.org/10.1111/biom.12139. The modification includes the FDR or FWER envelopes
#' (as specified by the argument \code{typeone} in \code{...}, passed to \code{\link{global_envelope_test}})
#' for the test of local correlations, i.e. multiple testing correction and
#' graphical illustration of the test results.
#'
#' Variograms are calculated using the package \pkg{geoR} and the local correlations using the R
#' package \pkg{locfit}. These packages should be installed to use \code{GET.localcor}.
#'
#' Currently the data is provided in the format of Vilodomat et al. (2014, Supporting information).
#' Additionally width and height of area represented by a data point can be provided, see the argument \code{data}.
#' This information is used for plotting purposes when plotting the output by \code{plot()}.
#'
#' Examples are provided in the vignette 'FDRenvelopes', see e.g. https://cran.r-project.org/package=GET.
#'
#' @param data A data.frame where the first two columns correspond to the values of the two random fields,
#' whose correlations are to be studied, and the third and fourth columns correspond to the x- and y-coordinates
#' where these random fields have been observed. In addition, the width and height of the pixels at each (x,y)
#' can be given in the fifth and sixth column. Warning: no checks for the data input.
#' @param Delta A smoothing parameter of the local correlation.
#' According to Vilodomat et al. (2014): Delta is a set of values for the proportion of neighbors
#' to consider for the smoothing step. No default. The user may have to experiment with different
#' values to find one suitable for their data.
#' @param nsim The number of resamples.
#' @param varying.bandwidth Logical, whether to use a varying bandwidth to calculate the local correlations or not.
#' See Vilodomat et al. (2014).
#' @param bandwidth.nn Nearest neighbor component of the smoothing parameter for varying bandwidth to
#' be passed to the argument \code{nn} of the function \code{lp} of the \pkg{locfit} package.
#' The user may have to experiment with different values to find one suitable for their data.
#' Default set to to 0.1 according to Vilodamat et al. (2014, supporting information).
#' @param bandwidth.h Non-varying bandwidth, to be passed to the argument \code{h} of
#' the function \code{lp} of \pkg{locfit}.
#' The user may have to experiment with different values to find one suitable for their data.
#' Default to 5.281 according to Vilodamat et al. (2014, supporting information).
#' @param maxk See \code{locfit} and \code{locfit.raw} of \pkg{locfit}. Default here to 300
#' following Vilodomat et al. (2014).
#' @param N_s If the number of observations is bigger than \code{N_s}, following Vilodomat et al. (2014)
#' a subsample of size N_s is taken every time when a variogram is calculated.
#' @param notest Logical. FALSE means that the test is done. TRUE allows to calculate only local
#' correlation for the data, which can be beneficial for choosing the bandwidths before running
#' the test. If TRUE, then only the observed local correlations will be returned.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' Note that for testing local correlations, it may often be preferable to use FDR control.
#' This can be specified by setting \code{typeone = "fdr"}, while the default is FWER control.
#' See \code{\link{global_envelope_test}} for defaults and available options.
#' @inheritParams graph.fanova
#' @inheritParams graph.flm
#' @importFrom stats cor
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @export
#'
#' @return A global envelope object (with possible additional classes), see description of main components 
#' in \code{\link{global_envelope}} (Value).
#' Additional attributes: \code{p_global} contains the Monte Carlo p-value for the global test of correlation.
#' \code{cor_global} and \code{cor_global_sim} contain the value of the correlation for data and permuted data,
#' respectively. If \code{savefuns = TRUE}, then \code{permutations} contain the permuted values of the first
#' random field according to Viladomat et al. (2014) procedure, and \code{cset} contains all the local correlations
#' for the data and permuted data in a \code{curve_set} object (see \code{\link{create_curve_set}}).
#'
#' @references
#' Viladomat, J., Mazumder, R., McInturff, A., McCauley, D.J. and Hastie, T. (2014). Assessing the significance of global and local correlations under spatial autocorrelation: A nonparametric approach. Biometrics 70, 409-418. doi: 10.1111/biom.12139
#'
#' Mrkvička, T., Myllymäki, M. (2023) False discovery rate envelopes. Statistics and Computing 33, 109. https://doi.org/10.1007/s11222-023-10275-7
#'
GET.localcor <- function(data, Delta, nsim = 1000, ...,
                         varying.bandwidth = FALSE,
                         bandwidth.nn = 0.1, bandwidth.h = 5.281, maxk = 300,
                         savefuns = FALSE, N_s = 1000,
                         mc.cores = 1L, mc.args = NULL, cl = NULL,
                         notest = FALSE) {
  if(nsim < 2) stop("Not a reasonable value of nsim.")
  if(!is.logical(varying.bandwidth)) stop("Invalid varying.bandwidth value.")
  if(!is.logical(savefuns)) stop("Invalid savefuns value.")
  if(!(is.numeric(Delta) & length(Delta)>0 & all(Delta>0))) stop("Invalid Delta.")
  # load the data X, Y and the coordinates of the N data locations:
  if(!(dim(data)[2] %in% c(4,6))) stop("data does not have four or six columns.")
  X <- data[,1]
  Y <- data[,2]
  lat <- data[,3]
  long <- data[,4]
  if(dim(data)[2] == 6) {
    d.x <- data[,5]
    d.y <- data[,6]
    widthheight <- TRUE
  }
  else {
    lat.x <- sort(unique(lat))
    long.y <- sort(unique(long))
    if(all(lat.x[-1] - lat.x[-length(lat.x)] == lat.x[2]-lat.x[1]) &
       all(long.y[-1] - long.y[-length(long.y)] == long.y[2]-long.y[1])) {
      widthheight <- TRUE
      d.x <- lat.x[2]-lat.x[1]
      d.y <- long.y[2]-long.y[1]
    }
    else widthheight <- FALSE
  }
  N <- length(X)

  # If N is too big, we take a subsample of size N_s every time we calculate a variogram:
  if(length(X) > N_s) {
    ids <- sample(N, N_s)
    X_s <- X[ids]
    long_s <- long[ids]
    lat_s <- lat[ids]
  } else {
    X_s <- X
    long_s <- long
    lat_s <- lat
    ids <- 1:N
  }

  dists <- dist(cbind(lat_s, long_s))

  # maximum distance for the variogram set at the 25% percentile of 
  # the distribution of pairs of distances:
  prctile <- quantile(dists, probs = 0.25)

  # variogram of variable X that will be used as target when doing the matching:
  target_variog <- geoR::variog(data = X_s, coords = cbind(long_s,lat_s), max.dist = prctile, option = "bin", messages = FALSE)

  # ALGORITHM:
  # It returns B random fields with the same autocorrelation 
  # as X but independent of Y, stored in permutations. The basis
  # to calculate B realizations of the null we are interested in.

  if(!notest) {
    loopfun1 <- function(i) {
      #print(c('i=', i))
      # random permutation of the values of X across locations:
      X.randomized <- sample(X, size = length(X), replace = FALSE)
      # smoothing and scaling step to match the target variogram:
      matching(X.randomized, long, lat, Delta, target_variog, prctile, ids, maxk)$hat.X.delta.star
    }
    if(is.null(cl)) permutations <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun1, mc.cores=mc.cores), mc.args))
    else permutations <- parLapply(cl, 1:nsim, loopfun1)
    permutations <- simplify2array(permutations)
  }

  # ASSESSING THE SINGLE PEARSON'S CORRELATION COEFFICIENT (GLOBAL CORRELATION):

  # observed global correlation:
  cor.global.obs <- as.vector(cor(X,Y))

  if(!notest) {
    # null distribution for the global correlation:
    cor.global <- cor(permutations,Y)

    # p-value:
    p.value.global <- estimate_p_value(abs(cor.global.obs), as.numeric(abs(cor.global)))
    # Same as 1 - sum(abs(cor.global) < abs(cor.global.obs)) / (nsim+1) = (1 + sum(abs(cor.global) >= abs(cor.global.obs))) / (nsim+1)
    #Vilodomat: p.value.global <- sum(abs(cor.global) > abs(cor.global.obs)) / nsim
  }

  # ASSESSING THE CORRELATION BETWEEN X AND Y IN A GIVEN NEIGHBORHOOD (LOCAL CORRELATIONS)

  # observed local correlations calculated using locfit:
  XY <- X*Y
  X2 <- X^2
  Y2 <- Y^2
  if(varying.bandwidth) {
    # nn: proportion of neighbors (varying bandwitdh). The default is 10% of N, but it may be changed by the user to find one more suitable for their data.
    fitX <- fitted(locfit::locfit(X ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
    fitY <- fitted(locfit::locfit(Y ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
    fitXY <- fitted(locfit::locfit(XY ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
    fitX2 <- fitted(locfit::locfit(X2 ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
    fitY2 <- fitted(locfit::locfit(Y2 ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
  } else {
    # h: non-varying bandwidth. The user may have to experiment with different values to find one suitable for their data.
    fitX <- fitted(locfit::locfit(X ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
    fitY <- fitted(locfit::locfit(Y ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
    fitXY <- fitted(locfit::locfit(XY ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
    fitX2 <- fitted(locfit::locfit(X2 ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
    fitY2 <- fitted(locfit::locfit(Y2 ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk)) 
  }
  loc.cor.obs <- (fitXY - fitX * fitY) / (sqrt((fitX2 - fitX ^ 2) * (fitY2 - fitY ^ 2)))

  if(!notest) {
    # loc.cor[i,] is the null distribution for the observed local correlation loc.cor.obs[i] at location i=1:N :
    loopfun2 <- function(k) {
      Xper <- permutations[,k]
      XYper <- Xper*Y
      Xper2 <- Xper^2
      if(varying.bandwidth) {
        fitXper <- fitted(locfit::locfit(Xper ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
        fitXYper <- fitted(locfit::locfit(XYper ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
        fitXper2 <- fitted(locfit::locfit(Xper2 ~ locfit::lp(long, lat, nn = bandwidth.nn, deg = 0), kern = "gauss", maxk = maxk))
      } else {
        fitXper <- fitted(locfit::locfit(Xper ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
        fitXYper <- fitted(locfit::locfit(XYper ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))
        fitXper2 <- fitted(locfit::locfit(Xper2 ~ locfit::lp(long, lat, h = bandwidth.h, deg = 0), kern = "gauss", maxk = maxk))   
      }
      # r_k
      (fitXYper - fitXper * fitY) / (sqrt((fitXper2 - fitXper ^ 2) * (fitY2 - fitY ^ 2)))
    }
    if(is.null(cl)) loc.cor <- do.call(mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args))
    else loc.cor <- parLapply(cl, 1:nsim, loopfun2)
    loc.cor <- simplify2array(loc.cor)

    if(widthheight)
      cset <- create_curve_set(list(r=data.frame(x=lat, y=long, width=d.x, height=d.y), obs=loc.cor.obs, sim_m=loc.cor))
    else
      cset <- create_curve_set(list(obs=loc.cor.obs, sim_m=loc.cor))

    res <- global_envelope_test(cset, ..., lower=-1, upper=1)

    attr(res, "p_global") <- p.value.global
    attr(res, "cor_global") <- cor.global.obs
    attr(res, "cor_global_sim") <- cor.global
    if(savefuns) {
      attr(res, "permutations") <- permutations
      attr(res, cset) <- cset
    }
  }
  else {
    if(widthheight)
      res <- create_curve_set(list(r=data.frame(x=lat, y=long, width=d.x, height=d.y), obs=loc.cor.obs))
    else
      res <- create_curve_set(list(r=seq_along(lat), obs=loc.cor.obs))
  }
  res
}
