#-----------------------------------------------------------#
# n sample test of correspondence of distribution functions #
#-----------------------------------------------------------#

# Ecdf means
# y is a vector for which groups gives grouping
#' @importFrom stats ecdf
ecdfmeans.m <- function(x, groups, r) {
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  # res <- sapply(ecdf.ls, FUN = function(x) { x(r) }, simplify=TRUE)
  # Names
  k <- nlevels(groups)
  gnames <- levels(groups)
  cont <- matrix(0, nrow=length(r), ncol=k)
  cont.names <- vector(length=k)
  for(i in 1:k){
    cont[,i] <- ecdf.ls[[i]](r) #res[,i] ## CHECK: Here rather  ecdf.ls[[i]](r) ??
    cont.names[i] <- paste("ECDF:", gnames[i])
  }
  colnames(cont) <- cont.names
  cont
}

# Ecdf contrasts
# y is a vector for which groups gives grouping
#' @importFrom stats ecdf
ecdfcontrasts.m <- function(x, groups, r) {
  k <- nlevels(groups)
  gnames <- levels(groups)
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  cont <- matrix(0, nrow=length(r), ncol=k*(k-1)/2)
  cont.names <- vector(length=k*(k-1)/2)
  counter <- 1
  for(i in 1:(k-1)) for(j in (i+1):k) {
    cont[, counter] <- ecdf.ls[[i]](r) - ecdf.ls[[j]](r)
    cont.names[counter] <- paste0("DIFF: ", paste(gnames[i], gnames[j], sep="-"))
    counter <- counter+1
  }
  colnames(cont) <- cont.names
  cont
}

# Density means
#' @importFrom stats density
denmeans.m <- function(x, groups, r, density.args) {
  k <- nlevels(groups)
  gnames <- levels(groups)
  cont <- matrix(0, nrow=length(r), ncol=k)
  cont.names <- vector(length=k)
  x.ls <- by(x, INDICES=groups,  FUN=function(x) x, simplify=FALSE)
  res <- sapply(x.ls,
              FUN = function(x) {
                do.call(density,
                        c(list(x = x,
                               from = min(r),
                               to = max(r),
                               n = length(r)),
                          density.args))$y
              },
              simplify=TRUE)
  for(i in 1:k) {
    cont[,i] <- res[,i]
    cont.names[i] <- paste0("DEN: ", gnames[i])
  }
  colnames(cont) <- cont.names
  cont
}

# Quantile regression
qreg.m <- function(x, groups, r, rq.args) {
  k <- nlevels(groups)
  gnames <- levels(groups)
  cont <- matrix(0, nrow=length(r), ncol = k-1)
  cont.names <- vector(length = k-1)
  fit <- do.call(quantreg::rq,
                 c(list(formula =  x ~ groups,
                        tau =r), rq.args)
  )
  cont <- matrix(coef( fit )[2:k,],
                 ncol = k - 1,
                 byrow = TRUE)
  cont.names <- vector(length = k-1)
  for(i in 2:k){
    cont.names[i-1] <- paste0("QR: ", gnames[i])
  }
  colnames(cont) <- cont.names
  cont
}

#' @importFrom stats ecdf
#' @importFrom stats approxfun
qq.m <- function(x, groups, r, approxfun.args, shift=FALSE) {
  k <- nlevels(groups)
  gnames <- levels(groups)
  cont <- matrix(0, nrow=length(r), ncol=k*(k-1)/2)
  cont.names <- vector(length=k*(k-1)/2)
  if(!shift) statname <- "QQ: "
  else statname <- "SHIFT: "
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  x.ls <- by(x, INDICES=groups, FUN=function(x) x, simplify=FALSE)
  counter <- 1
  for(i in 1:(k-1)) for(j in (i+1):k) {
    # A function that approximates quantiles
    quant <- do.call(approxfun,
                    c(list(x = seq(0, 1, length = length(x.ls[[j]])),
                           y = sort(x.ls[[j]]),
                           yleft = NA,
                           yright = NA), approxfun.args)
    )
    if(!shift) cont[, counter] <- quant(ecdf.ls[[i]](r))
    else cont[, counter] <- quant(ecdf.ls[[i]](r)) - r
    cont.names[counter] <- paste0(statname, gnames[i], "-", gnames[j])
    counter <- counter+1
  }
  colnames(cont) <- cont.names
  cont
}


#' Graphical n sample test of correspondence of distribution functions
#'
#' Compare the distributions of two (or more) samples.
#'
#'
#' A global envelope test can be performed to investigate whether the n distribution functions
#' differ from each other and how do they differ. This test is a generalization of
#' the two-sample Kolmogorov-Smirnov test with a graphical interpretation.
#' We assume that the observations in the sample \eqn{i}{i} are an i.i.d. sample from the distribution
#' \eqn{F_i(r), i=1, \dots, n,}{F_i(r), i=1, ..., n,}
#' and we want to test the hypothesis
#' \deqn{F_1(r)= \dots = F_n(r).}{F_1(r)= ... = F_n(r).}
#' If \code{contrasts = FALSE} (default), then the default test statistic ("ECDF") is taken to be
#' \deqn{\mathbf{T} = (\hat{F}_1(r), \dots, \hat{F}_n(r))}{T = (\hat{F}_1(r), \dots, \hat{F}_n(r))}
#' where \eqn{\hat{F}_i(r) = (\hat{F}_i(r_1), \dots, \hat{F}_i(r_k))}{\hat{F}_i(r) = (\hat{F}_i(r_1), ..., \hat{F}_i(r_k))}
#' is the ecdf of the \eqn{i}{i}th sample evaluated at argument values
#' \eqn{r = (r_1,\dots,r_k)}{r = (r_1, ...,r_k)}.
#'
#' Another possibility is given by \code{stat = "DIFF"}, and then the test statistic is
#' still based on the ECDFs and constructed from all pairwise differences,
#' \deqn{\mathbf{T} = (\hat{F}_1(r)-\hat{F}_2(r), \hat{F}_1(r)-\hat{F}_3(r), \dots, \hat{F}_{n-1}(r)-\hat{F}_n(r))}{T = (\hat{F}_1(r)-\hat{F}_2(r), ..., \hat{F}_{n-1}(r)-\hat{F}_n(r))}
#' The choices \code{contrasts = TRUE} and \code{stat = "ECDF"} lead to the same test statistic.
#' For other (or multiple) values of \code{stat}, the argument \code{contrasts} is ignored.
#'
#' All the options as the test statistics are the following:
#' \enumerate{
#' \item \code{"ECDF"}: The ECDFs of the n-samples, as specified above
#' \item \code{"DIFF"}: The pairwise differences between the ECDFs, as specified above
#' \item \code{"DEN"}: The kernel estimated density functions of the n-samples as the test statistic
#' \item \code{"QQ"}: The pairwise comparisons of empirical quantiles
#' \item \code{"SHIFT"} The de-trended QQ-plot (shift plot)
#' \item \code{"QR"}: The quantile regression process, i.e. the \eqn{\beta}{beta}-coefficients of the quantile regression.
#' By default, the reference category of this test statistic is the first sample.
#' }
#' The test statistics are described in detail in Konstantinou et al. (2024).
#'
#' The simulations under the null hypothesis that the distributions are the same are obtained
#' by permuting the individuals of the groups. The default number of permutation, if \code{nsim} is not specified,
#' is \eqn{n \cdot 1000-1}{n*1000-1} for the case \code{contrasts = FALSE} and
#' \eqn{(n \cdot (n-1)/2) \cdot 1000 - 1}{(n*(n-1)/2)*1000 - 1} for the case \code{contrasts = TRUE},
#' where \eqn{n}{n} is the length of \eqn{x}{x}.
#'
#' @inheritParams graph.fanova
#' @inheritParams graph.flm
#' @param x A list of numeric vectors, one for each sample.
#' @param r The sequence of argument values at which the test functions are to be compared.
#' The default is 100 equally spaced values between the minimum and maximum over all groups.
#' @param tau The sequence of argument values for the QR test statistic.
#' The default values are 100 equally spaced values between 0.1 and 0.9.
#' @param stat Character string indicating which test statistic to be used. See details.
#' @param density.args  A named list of additional arguments to be passed for the estimation of the test statistic \code{"DEN"}.
#' For more details see \code{\link[stats]{density}}.
#' @param approxfun.args  A named list of additional arguments to be passed for the estimation of the the test statistic \code{"QQ"}.
#' For more details see \code{\link[stats]{approxfun}}.
#' @param rq.args A named list of additional arguments to be passed for the estimation of the test statistic \code{"QR"}.
#' For more details see the function \code{rq} of \pkg{quantreq}.
#' @export
#' @references Konstantinou K., Mrkvička T. and Myllymäki M. (2024) Graphical n-sample tests of correspondence of distributions. arXiv:2403.01838 [stat.ME] https://doi.org/10.48550/arXiv.2403.01838
#' @aliases GET.necdf
#' @examples
#' if(require("fda", quietly=TRUE)) {
#'   # Heights of boys and girls at age 10
#'   f.a <- growth$hgtf["10",] # girls at age 10
#'   m.a <- growth$hgtm["10",] # boys at age 10
#'   # Empirical cumulative distribution functions
#'   plot(ecdf(f.a))
#'   plot(ecdf(m.a), col='grey70', add=TRUE)
#'   # Create a list of the data
#'   fm.list <- list(Girls=f.a, Boys=m.a)
#'   \donttest{
#'   res <- GET.distrequal(fm.list)
#'   plot(res)
#'   # If you want to change the labels:
#'   plot(res, scales = "free", labels = c("Girls", "Boys"))
#'   # If you want to change the x-label (y-label similarly):
#'   require("ggplot2")
#'   myxlab <- substitute(paste(italic(i), " (", j, ")", sep = ""),
#'                        list(i = "x", j = "Height in cm"))
#'   plot(res, scales = "free") + xlab(myxlab)
#'   # Use instead the test statistics QQ and DEN
#'   res <- GET.distrequal(fm.list, stat = c("QQ", "DEN"))
#'   plot(res, scales = "free")
#'   }
#'   \dontshow{
#'   # The test with a lower number of simulations
#'   res <- GET.distrequal(fm.list, nsim=4, alpha=0.2)
#'   plot(res)
#'   res <- GET.distrequal(fm.list, stat = c("QQ", "DEN"), nsim=4, alpha=0.2)
#'   plot(res, scales = "free")
#'   }
#'
#'   # Heights of boys and girls at age 14
#'   f.a <- growth$hgtf["14",] # girls at age 14
#'   m.a <- growth$hgtm["14",] # boys at age 14
#'   # Empirical cumulative distribution functions
#'   plot(ecdf(f.a))
#'   plot(ecdf(m.a), col='grey70', add=TRUE)
#'   # Create a list of the data
#'   fm.list <- list(Girls=f.a, Boys=m.a)
#'   \donttest{
#'   res <- GET.distrequal(fm.list)
#'   plot(res) + xlab(myxlab)
#'   res <- GET.distrequal(fm.list, stat = c("QQ", "DEN"))
#'   plot(res, scales = "free") + xlab(myxlab)
#'   }
#'   \dontshow{
#'   # The test with a lower number of simulations
#'   res_m <- GET.distrequal(fm.list, nsim=4, alpha=0.2)
#'   plot(res_m)
#'   res_c <- GET.distrequal(fm.list, stat = c("QQ", "DEN"), nsim=4, alpha=0.2)
#'   plot(res_c, scales = "free")
#'   }
#' }
#' if(require("datasets", quietly=TRUE)) {
#'   data("iris")
#'   virginica <- subset(iris, Species == "virginica")
#'   setosa <- subset(iris, Species == "setosa")
#'   versicolor <- subset(iris, Species == "versicolor")
#'   \donttest{
#'   res <- GET.distrequal(x = list(virginica = virginica$Sepal.Length,
#'                                  setosa = setosa$Sepal.Length,
#'                                  versicolor = versicolor$Sepal.Length),
#'                         stat =  c("QQ", "DEN"))
#'   plot(res, scales = "free", ncol = 3)
#'   }
#'   \dontshow{
#'   res <- GET.distrequal(x = list(virginica = virginica$Sepal.Length,
#'                                  setosa = setosa$Sepal.Length,
#'                                  versicolor = versicolor$Sepal.Length),
#'                         stat =  c("QQ", "DEN"), nsim=4, alpha=0.2)
#'   plot(res, scales = "free", ncol = 3)
#'   }
#' }
GET.distrequal <- function(x, stat = "ECDF", nsim,
                           r = seq(min(unlist((lapply(x, min)))), max(unlist((lapply(x, max)))), length=100),
                           tau = seq(0.1, 0.9, length=100),
                           contrasts = FALSE,
                           GET.args = NULL,
                           density.args = NULL,
                           approxfun.args = NULL,
                           rq.args=NULL,
                           savefuns = FALSE, ...) {
  if(!is.list(x) && length(x)<2) stop("At least two groups should be provided.")
  x.lengths <- as.numeric(lapply(x, FUN = length))
  if(!is.null(names(x))) groups <- rep(names(x), times=x.lengths)
  else groups <- rep(seq_along(x), times=x.lengths)
  groups <- factor(groups, levels=unique(groups))
  gnames <- levels(groups)
  k <- nlevels(groups)
  allowed_test.stat <- c("ECDF", "DEN", "DIFF", "QQ", "SHIFT", "QR")
  stat <- match.arg(toupper(stat),
                    choices = allowed_test.stat,
                    several.ok = TRUE)
  if(length(stat)>1) contrasts <- FALSE
  GET.args <- list(...)
  # The default nstep
  if(!("nstep" %in% GET.args)) {
    if(length(stat)==1) nstep <- 1
    else nstep <- 2
    GET.args <- c(list(nstep=nstep), GET.args)
  }
  if("QR" %in% stat && (min(tau)<0 || max(tau)>1))
    stop("Invalid discretization of the quantile regression stat! Please provide tau values that are >=0 and <=1.")

  # Setting the 'fun'
  choose_fun <- function(stat_type) {
    switch(stat_type,
                ECDF = {
                  if(!contrasts) ecdfmeans.m
                  else ecdfcontrasts.m
                },
                DIFF = { ecdfcontrasts.m },
                DEN  = { function(x, groups, r) denmeans.m(x, groups, r, density.args) },
                QR   = { function(x, groups, r) qreg.m(x, groups, tau, rq.args) },
                QQ   = { function(x, groups, r) {
                  qq.m(x, groups, r, approxfun.args, shift=FALSE)
                }},
                SHIFT = { function(x, groups, r) {
                  qq.m(x, groups, r, approxfun.args, shift=TRUE)
                }})
  }

  x <- unlist(x)

  # Observed difference between the ecdfs, and other test statistics
  obs <- vector(mode="list", length=length(stat))
  for(i in seq_along(stat)) { # Go through all chosen test statistics
    fun <- choose_fun(stat[i])
    obs[[i]] <- fun(x, groups, r)
  }

  # How many csets will we have?
  l <- sum(sapply(obs, FUN = ncol, simplify=TRUE))
  # How many simulations?
  if(missing(nsim)) {
    nsim <- l * 1000 - 1
    message("Creating ", nsim, " permutations.\n", sep="")
  }

  # Simulations by permuting to which groups each value belongs to
  sim <- replicate(nsim,
                   expr = {
                     #fun(x, sample(groups, size=length(groups), replace=FALSE), r), simplify = "array")
                     permgroups <- sample(groups, size=length(groups), replace=FALSE)
                     tmp <- vector(mode="list", length=length(stat))
                     for(i in seq_along(stat)) { # Go through all chosen test statistics
                       fun <- choose_fun(stat[i])
                       tmp[[i]] <- fun(x, permgroups, r)
                     }
                     tmp
                   }, simplify = FALSE)
  complabels <- unlist(lapply(obs, FUN = colnames))

  csets <- vector("list", length=l)
  counter <- 1
  for(j in seq_along(obs)) {
    for(i in 1:ncol(obs[[j]])) {
      csets[[counter]] <- create_curve_set(
        list(r = r, # OR tau
             obs = obs[[j]][,i],
             sim_m = sapply(sim, FUN = function(x) x[[j]][,i], simplify=TRUE)
             ))
      counter <- counter + 1
    }
  }
  names(csets) <- complabels

  # FWER or FDR envelope
  res <- do.call(global_envelope_test, c(list(curve_sets=csets, alternative="two.sided"), GET.args))
  if(length(stat) == 1 && stat == "ECDF") {
    if(!contrasts)
      res <- envelope_set_labs(res, xlab = expression(italic(x)),
                               ylab = expression(italic(hat(F)(x))))
    else
      res <- envelope_set_labs(res, xlab = expression(italic(x)),
                               ylab = expression(italic(hat(F)[i](x)-hat(F)[j](x))))
  }
  attr(res, "contrasts") <- contrasts
  attr(res, "labels") <- complabels
  attr(res, "method") <- paste("Graphical n-sample test" )
  if(savefuns) attr(res, "simfuns") <- csets
  attr(res, "call") <- match.call()
  res
}
