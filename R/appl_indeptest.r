#' Test of independence based on cumulative distribution function
#'
#' Permutation-based test of independence in a bivariate vector
#' using the empirical joint cumulative distribution function as
#' the test statistic.
#'
#'
#' Permutation-based test of independence in a bivariate sample,
#' based on empirical joint cumulative distribution function
#' computed on a grid of \code{ngrid[1]} times \code{ngrid[2]} arguments.
#' The grid points are chosen according to the quantiles of the
#' marginal distributions.
#'
#' If the observed data are the pairs \eqn{\{(X_1, Y_1), \ldots, (X_n, Y_n)\}}{{(X_1, Y_1), ..., (X_n, Y_n)}},
#' the permutations are obtained by randomly permuting the values
#' in the second marginal, i.e. \eqn{\{(X_1, Y_{\pi(1)}), \ldots, (X_n, Y_{\pi(n)})\}}{{(X_1, Y_{pi(1)}), ..., (X_n, Y_{pi(n)})}}.
#'
#' The test itself is performed using the global envelope test of the chosen version,
#' see the argument \code{type} of \code{\link{global_envelope_test}}.
#'
#' @param X A matrix with n rows and 2 columns. Each row contains
#' one bivariate observation.
#' @param ngrid Vector with two elements, giving the number of grid
#' points to be used in the test statistic for each of the two marginals.
#' The default is 20 in each marginal.
#' @param nsim The number of random permutations used.
#' @param seq.x For the first marginal, the values at which the
#' empirical cumulative distribution function will be evaluated.
#' If NULL (the default), sequence of quantiles will be used,
#' equidistant in terms of probability.
#' @param seq.y For the second marginal, the values at which the
#' empirical cumulative distribution function will be evaluated.
#' If NULL (the default), sequence of quantiles will be used,
#' equidistant in terms of probability.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' In particularly, \code{alpha} specifies the nominal significance level of the test,
#' and \code{type} the type of the global envelope test.
#' @export
#' @references Dvořák, J. and Mrkvička, T. (2022). Graphical tests of independence for general distributions. Computational Statistics 37, 671--699.
#' @examples
#' # Generate sample data
#' data <- matrix(rnorm(n=200), ncol=2) %*% matrix(c(1,1,0,1), ncol=2)
#' plot(data)
#'
#' # Compute the CDF test and plot the significant regions
#' \donttest{res <- GET.cdf(data, ngrid=c(20,15), nsim=1999)}
#' \dontshow{res <- GET.cdf(data, ngrid=c(20,15), nsim=19)}
#'
#' # Graphical output of GET
#' plot(res)
#' plot(res, what=c("lo.sign", "hi.sign"))
#'
#' # Extract the p-value
#' attr(res,"p")
GET.cdf <- function(X, ngrid=c(20,20), nsim=999, seq.x=NULL, seq.y=NULL, ...) {

  X <- as.matrix(X)
  if(ncol(X) != 2) stop("Invalid X. X should have two columns.")
  # Number of observations
  n <- nrow(X)

  # Choosing the bins if seq.x or seq.y is not provided
  if(is.null(seq.x)) { seq.x <- as.numeric(quantile(X[,1], probs=seq(from=1/ngrid[1],to=1, length.out=ngrid[1]))) }
  else if(!is.numeric(seq.x) | !is.vector(seq.x)) stop("seq.x should be a numeric vector.")
  if(is.null(seq.y)) { seq.y <- as.numeric(quantile(X[,2], probs=seq(from=1/ngrid[2],to=1, length.out=ngrid[2]))) }
  else if(!is.numeric(seq.y) | !is.vector(seq.y)) stop("seq.y should be a numeric vector.")

  # Computing cdf values in observed data
  F.obs <- matrix(NA, ncol=length(seq.y),nrow=length(seq.x))
  for(k in 1:length(seq.x)) {
    for(j in 1:length(seq.y)) {
      F.obs[k,j] <- sum( (X[,1] <= seq.x[k]) & (X[,2] <= seq.y[j]) ) / n
    }
  }

  # Computing cdf values in permuted data
  F.perm <- array(NA, dim=c(length(seq.x),length(seq.y),nsim))
  for (i in 1:nsim){
    X.perm <- cbind(X[,1],X[sample.int(n),2])
    for (k in 1:length(seq.x)){
      for (j in 1:length(seq.y)){
        F.perm[k,j,i] <- sum( (X.perm[,1] <= seq.x[k]) & (X.perm[,2] <= seq.y[j]) ) / n
      }
    }
  }

  # Performing 2D global envelope test
  im.set <- create_image_set(list(r=list(1:length(seq.x),1:length(seq.y)), obs=F.obs, sim_m=F.perm))
  erl2d <- global_envelope_test(im.set, ...)
  attr(erl2d, "method") <- "cdf-based permutation test"

  # Returning the output of the global envelope test
  return(erl2d)
}

#' Test of independence in a 2D contingency table
#'
#' Permutation-based test of independence in a 2D contingency table,
#' using the matrix of observed counts as the test statistic.
#'
#'
#' Permutation-based test of independence in a 2D contingency table,
#' using the matrix of observed counts as the test statistic.
#'
#' If the observed data are the pairs \eqn{\{(X_1, Y_1), \ldots, (X_n, Y_n)\}}{{(X_1, Y_1), ..., (X_n, Y_n)}},
#' the permutations are obtained by randomly permuting the values
#' in the second marginal, i.e. \eqn{\{(X_1, Y_{\pi(1)}), \ldots, (X_n, Y_{\pi(n)})\}}{{(X_1, Y_{pi(1)}), ..., (X_n, Y_{pi(n)})}}.
#'
#' The test itself is performed using the global envelope test in the chosen version.
#' Text output can be printed in the console by typing the object name.
#' The cells in which the observed value exceeds the upper envelope printed in red,
#' and cells in which the observed value is lower than the lower
#' envelope printed in cyan. Standard output of the global envelope
#' test is also returned and can be plotted or analyzed accordingly.
#' @inheritParams GET.cdf
#' @export
#' @references Dvořák, J. and Mrkvička, T. (2022). Graphical tests of independence for general distributions. Computational Statistics 37, 671--699.
#' @examples
#' # Generate sample data:
#' data <- matrix(c(sample(4, size=100, replace=TRUE), sample(2, size=100, replace=TRUE)), ncol=2)
#' data[,2] <- data[,2] + data[,1]
#'
#' # Observed contingency table (with row names and column names)
#' table(data[,1], data[,2])
#'
#' # Permutation-based envelope test
#' \donttest{res <- GET.contingency(data, nsim=999)}
#' \dontshow{res <- GET.contingency(data, nsim=19)}
#' res
#' plot(res)
#'
#' # Extract the p-value
#' attr(res,"p")
GET.contingency <- function(X, nsim=999, ...){

  X <- as.matrix(X)
  if(ncol(X) != 2) stop("Invalid X. X should have two columns.")
  # Number of observations
  n <- nrow(X)

  # Number of categories in each marginal
  k1 <- length(unique(X[,1]))
  k2 <- length(unique(X[,2]))

  # Contingency table for observed data
  tab <- table(X[,1],X[,2])

  # Contingency tables for permuted data
  sim_m <- array(NA, dim=c(k1,k2,nsim))
  for (i in 1:nsim){
    sim_m[,,i] <- table(X[,1],X[sample(n),2])
  }

  # Performing 2D global envelope test
  cs <- create_image_set(list(obs=tab, sim_m=sim_m))
  erl2d <- global_envelope_test(cs, ...)

  # Rearranging the output of the global envelope test to match the output of "table" (by print)
  aux <- erl2d$x
  erl2d$x <- erl2d$y
  erl2d$y <- aux
  aux <- matrix(erl2d$obs,nrow=k1,byrow=F)[k1:1,]
  erl2d$obs <- as.vector(aux)
  aux <- matrix(erl2d$central,nrow=k1,byrow=F)[k1:1,]
  erl2d$central <- as.vector(aux)
  aux <- matrix(erl2d$hi,nrow=k1,byrow=F)[k1:1,]
  erl2d$hi <- as.vector(aux)
  aux <- matrix(erl2d$lo,nrow=k1,byrow=F)[k1:1,]
  erl2d$lo <- as.vector(aux)

  # Returning the output of the global envelope test
  class(erl2d) <- c("GET_contingency", class(erl2d))
  attr(erl2d, "ncategories") <- c(k1, k2)
  attr(erl2d, "method") <- "Permutation test for contingency tables"
  erl2d
}

#' Print method for the class 'GET_contingency'
#'
#' @param x A 'GET_contingency' object
#' @param ... Ignored.
#' @export
print.GET_contingency <- function(x, ...) {
  k1 <- attr(x, "ncategories")[1]
  # Visualization of significant cells
  cat("Observed counts with significant cells highlighted in color:\n")
  obs <- matrix(x$obs,nrow=k1)[k1:1,]
  hi <- matrix(x$hi,nrow=k1)[k1:1,]
  lo <- matrix(x$lo,nrow=k1)[k1:1,]
  len <- 1 + nchar(max(max(obs),max(hi),max(lo)))
  for(j in 1:length(obs[,1])) {
    a <- ""
    for(i in 1:length(obs[1,])) {
      l <- nchar(obs[j,i])
      if(l < len){
        for(k in 1:(len-l)){
          a <- paste(a,"")
        }
      }
      if(obs[j,i] < lo[j,i]) { a <- paste(a,crayon::cyan(obs[j,i])) }
      if(obs[j,i] > hi[j,i]) { a <- paste(a,crayon::red(obs[j,i])) }
      if((obs[j,i] >= lo[j,i]) & (obs[j,i] <= hi[j,i])) { a <- paste(a,obs[j,i]) }
    }
    cat(a)
    cat("\n")
  }
  NextMethod()
}

#----------------------------------------------------------------------------#
# Test of independence in a bivariate sample based on quantile-quantile plot #
#----------------------------------------------------------------------------#

# Estimate the intensity function in the Q-Q representation
#' @importFrom stats dnorm pnorm
qq.estimate <- function(X, ngrid=c(64,64), sigma=NULL, atoms.x=NULL, atoms.y=NULL) {

  X <- as.matrix(X)
  if(ncol(X) != 2) stop("Invalid X. X should have two columns.")
  # Number of observations
  n <- nrow(X)

  # Default bw
  if(is.null(sigma)) { sigma <- n^(-1/6)*sqrt(1/12) }

  # Creating a point pattern from the Q-Q representation of the data
  X.pp <- suppressWarnings(spatstat.geom::ppp(x=rank(X[,1])/n, y=rank(X[,2])/n))

  # Two-dimensional estimates of the intensity function - observed data
  X.dens <- spatstat.core::density.ppp(X.pp, dimyx=c(ngrid[2],ngrid[1]), sigma=sigma)

  # One-dimensional estimates - observed data
  if(!is.null(atoms.y)) {
    wh.y <- rep(NA, times=length(atoms.y)) # which row is central
    hm.y <- rep(NA, times=length(atoms.y)) # how many rows
    for(i in 1:length(atoms.y)) {
      wh.y[i] <- X.pp$y[which(X[,2]==atoms.y[i])[1]]*ngrid[2]
      hm.y[i] <- sum(X[,2]==atoms.y[i])/n * ngrid[2]/2
    }
    grid.y <- X.dens$yrow
    bw <- attr(X.dens,"sigma")
    edgewt.y <- pnorm((1-grid.y)/bw) - pnorm((0-grid.y)/bw)
    for(i in 1:length(atoms.y)) {
      y.val <- atoms.y[i]
      ind <- X[,2] == y.val
      pts <- X.pp$x[ind]
      repl <- rep(NA, times=length(grid.y))
      for(j in 1:length(grid.y)) {
        t1 <- grid.y[j]
        weights1 <- dnorm(x=pts, mean=t1, sd=bw)
        repl[j] <- sum(weights1)
      }
      repl <- (repl/edgewt.y)/(sum(ind)/n) # scaling assuring the image values integrate to the number of observations
      for(k in round(wh.y[i]-hm.y[i]+1):round(wh.y[i]+hm.y[i])) {
        X.dens$v[k,] <- repl
      }
    }
  }
  if(!is.null(atoms.x)) {
    wh.x <- rep(NA, times=length(atoms.x)) # which column is central
    hm.x <- rep(NA, times=length(atoms.x)) # how many columns
    for(i in 1:length(atoms.x)) {
      wh.x[i] <- X.pp$x[which(X[,1]==atoms.x[i])[1]]*ngrid[1]
      hm.x[i] <- sum(X[,1]==atoms.x[i])/n * ngrid[1]/2
    }
    grid.x <- X.dens$xcol
    bw <- attr(X.dens,"sigma")
    edgewt.x <- pnorm((1-grid.x)/bw) - pnorm((0-grid.x)/bw)
    for(i in 1:length(atoms.x)) {
      x.val <- atoms.x[i]
      ind <- X[,1] == x.val
      pts <- X.pp$y[ind]
      repl <- rep(NA, times=length(grid.x))
      for(j in 1:length(grid.x)) {
        t1 <- grid.x[j]
        weights1 <- dnorm(x=pts, mean=t1, sd=bw)
        repl[j] <- sum(weights1)
      }
      repl <- (repl/edgewt.x)/(sum(ind)/n) # scaling assuring the image values integrate to the number of observations
      for(k in round(wh.x[i]-hm.x[i]+1):round(wh.x[i]+hm.x[i])) {
        X.dens$v[,k] <- repl
      }
    }
  }

  # Zero-dimensional estimates - observed data
  if(!is.null(atoms.x) & !is.null(atoms.y)) {
    for(i in 1:length(atoms.x)) {
      for(j in 1:length(atoms.y)) {
        x.val <- atoms.x[i]
        y.val <- atoms.y[j]
        npts <- sum((X[,1] == x.val) & (X[,2] == y.val))
        area <- 2*hm.x[i]*2*hm.y[j]/(ngrid[1]*ngrid[2])
        int <- npts/area
        X.dens[round(wh.y[j]-hm.y[j]+1):round(wh.y[j]+hm.y[j]),round(wh.x[i]-hm.x[i]+1):round(wh.x[i]+hm.x[i])] <- int
      }
    }
  }

  return(X.dens)
}

#' Test of independence based on the smoothed Q-Q plot
#'
#' Permutation-based test of independence in a bivariate vector
#' using the smoothed Q-Q plot as the test statistic.
#'
#'
#' Permutation-based test of independence in a bivariate sample,
#' based on Q-Q representation and estimate of the intensity function
#' computed on a regular grid of \code{ngrid[1]} times \code{ngrid[2]} points.
#'
#' If the observed data are the pairs \eqn{\{(X_1, Y_1), \ldots, (X_n, Y_n)\}}{{(X_1, Y_1), ..., (X_n, Y_n)}},
#' the permutations are obtained by randomly permuting the values
#' in the second marginal, i.e. \eqn{\{(X_1, Y_{\pi(1)}), \ldots, (X_n, Y_{\pi(n)})\}}{{(X_1, Y_{pi(1)}), ..., (X_n, Y_{pi(n)})}}.
#'
#' The test itself is performed using the global envelope test in the chosen version.
#'
#' @inheritParams GET.cdf
#' @param ngrid Vector with two elements, giving the number of grid
#' points to be used in the test statistic for each of the two marginals.
#' The default is 64 in each marginal.
#' @param sigma Standard deviation of the smoothing kernel to be 
#' used for smoothing the Q-Q plot when computing the test statistic.
#' If NULL, sensible default value is used based on the number of observations.
#' @param atoms.x Vector specifying atomic values in the first marginal. See Examples.
#' @param atoms.y Vector specifying atomic values in the second marginal. See Examples.
#' @export
#' @references Dvořák, J. and Mrkvička, T. (2022). Graphical tests of independence for general distributions. Computational Statistics 37, 671--699.
#' @examples
#' # Generate sample data
#' data <- matrix(rnorm(n=200), ncol=2) %*% matrix(c(1,1,0,1), ncol=2)
#' plot(data)
#'
#' # Compute the QQ test and plot the significant regions
#' \donttest{res <- GET.qq(data, ngrid=c(30,20), nsim=999)}
#' \dontshow{res <- GET.qq(data, ngrid=c(30,20), nsim=19)}
#' plot(res)
#' # Original graphical output of GET
#' plot(res)
#' # Extract the p-value
#' attr(res,"p")
#'
#' # With atoms, independent
#' data <- cbind(rnorm(n=100), sample(4, size=100, replace=TRUE))
#' plot(data)
#' \donttest{res <- GET.qq(data, nsim=999, atoms.y=c(1,2,3,4))}
#' \dontshow{res <- GET.qq(data, nsim=19, atoms.y=c(1,2,3,4))}
#' plot(res)
#' # Extract the p-value
#' attr(res,"p")
#'
#' # With atoms, dependent
#' data <- cbind(sort(rnorm(n=100)), sort(sample(4, size=100, replace=TRUE)))
#' plot(data)
#' \donttest{res <- GET.qq(data, nsim=999, atoms.y=c(1,2,3,4))}
#' \dontshow{res <- GET.qq(data, nsim=19, atoms.y=c(1,2,3,4))}
#' plot(res)
#' # Original graphical output of GET
#' plot(res)
#' # Extract the p-value
#' attr(res,"p")
#'
#' # Atoms in both variables
#' data <- cbind(rnorm(n=100), rnorm(n=100)) %*% matrix(c(1,1,0,1), ncol=2)
#' data[,1][data[,1]<=-1] <- -1
#' data[,2][data[,2]<=-0.5] <- -0.5
#' plot(data)
#' 
#' # Perform the test
#' \donttest{es <- GET.qq(data, nsim=999, atoms.x=c(-1), atoms.y=c(-0.5), sigma=NULL)}
#' \dontshow{es <- GET.qq(data, nsim=19, atoms.x=c(-1), atoms.y=c(-0.5), sigma=NULL)}
#' plot(res)
#' # Extract the p-value
#' attr(res,"p")
GET.qq <- function(X, ngrid=c(64,64), nsim=999, sigma=NULL, atoms.x=NULL, atoms.y=NULL, ...) {

  X <- as.matrix(X)
  if(ncol(X) != 2) stop("Invalid X. X should have two columns.")
  # Number of observations
  n <- nrow(X)

  # Creating a point pattern from the Q-Q representation of the data
  X.pp <- suppressWarnings(spatstat.geom::ppp(x=rank(X[,1])/n, y=rank(X[,2])/n))

  # Two-dimensional estimates of the intensity function - observed data
  X.dens <- qq.estimate(X, ngrid=ngrid, atoms.x=atoms.x, atoms.y=atoms.y, sigma=sigma)

  # Estimating the intensity function - permuted data
  sim_m <- array(NA, dim=c(ngrid[1],ngrid[2],nsim))
  for(j in 1:nsim) {
    X.perm <- cbind(X[,1],X[sample(n),2])
    X.dens.perm <- qq.estimate(X.perm, ngrid=ngrid, atoms.x=atoms.x, atoms.y=atoms.y, sigma=sigma)
    sim_m[,,j] <- t(X.dens.perm$v)
  }

  # Performing 2D global envelope test
  im.set <- create_image_set(list(obs=t(X.dens$v), sim_m=sim_m))
  erl2d <- global_envelope_test(im.set, ...)

  # Returning the output of the global envelope test and other features for plotting
  attr(erl2d, "method") <- "QQ-test of independence"
  return(erl2d)
}
