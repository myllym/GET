#-----------------------------------------------------------#
# n sample test of correspondence of distribution functions #
#-----------------------------------------------------------#

# Ecdf means
# y is a vector for which groups gives grouping
#' @importFrom stats ecdf
ecdfmeans.m <- function(x, groups, r) {
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  sapply(ecdf.ls, FUN = function(x) { x(r) }, simplify=TRUE)
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
    cont.names[counter] <- paste(gnames[i], gnames[j], sep="-")
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
#' If \code{contrasts = FALSE} (default), then the test statistic is taken to be
#' \deqn{\mathbf{T} = (\hat{F}_1(r), \dots, \hat{F}_n(r))}{T = (\hat{F}_1(r), \dots, \hat{F}_n(r))}
#' where \eqn{\hat{F}_i(r) = (\hat{F}_i(r_1), \dots, \hat{F}_i(r_k))}{\hat{F}_i(r) = (\hat{F}_i(r_1), ..., \hat{F}_i(r_k))}
#' is the ecdf of the \eqn{i}{i}th sample evaluated at argument values
#' \eqn{r = (r_1,\dots,r_k)}{r = (r_1, ...,r_k)}.
#' This is our recommended test function for the test.
#' Another possibility is given by \code{contrasts = TRUE}, and then the test statistic is contructed from
#' all pairwise differences,
#' \deqn{\mathbf{T} = (\hat{F}_1(r)-\hat{F}_2(r), \hat{F}_1(r)-\hat{F}_3(r), \dots, \hat{F}_{n-1}(r)-\hat{F}_n(r))}{T = (\hat{F}_1(r)-\hat{F}_2(r), ..., \hat{F}_{n-1}(r)-\hat{F}_n(r))}
#'
#' The simulations under the null hypothesis that the distributions are the same are obtained
#' by permuting the individuals of the groups. The default number of permutation, if nsim is not specified,
#' is n*1000 - 1 for the case \code{contrasts = FALSE} and
#' (n*(n-1)/2)*1000 - 1 for the case \code{contrasts = TRUE},
#' where n is the length of x.
#'
#' @param x A list of numeric vectors, one for each sample.
#' @param r The sequence of argument values at which the distribution functions are to be compared.
#' The default is 100 equally spaced values between the minimum and maximum over all groups.
#' @inheritParams graph.fanova
#' @export
#' @examples
#' if(require(fda, quietly=TRUE)) {
#'   # Heights of boys and girls at age 10
#'   f.a <- growth$hgtf["10",] # girls at age 10
#'   m.a <- growth$hgtm["10",] # boys at age 10
#'   # Empirical cumulative distribution functions
#'   plot(ecdf(f.a))
#'   plot(ecdf(m.a), col='grey70', add=TRUE)
#'   # Create a list of the data
#'   fm.list <- list(Girls=f.a, Boys=m.a)
#'   \donttest{
#'   res_m <- GET.necdf(fm.list)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE)
#'   plot(res_c)
#'   }
#'   \dontshow{
#'   # The test with lower number of simulations
#'   res_m <- GET.necdf(fm.list, nsim=4, alpha=0.2)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE, nsim=4, alpha=0.2)
#'   plot(res_c)
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
#'   res_m <- GET.necdf(fm.list)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE)
#'   plot(res_c)
#'   }
#'   \dontshow{
#'   # The test with lower number of simulations
#'   res_m <- GET.necdf(fm.list, nsim=4, alpha=0.2)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE, nsim=4, alpha=0.2)
#'   plot(res_c)
#'   }
#' }
GET.necdf <- function(x, r = seq(min(unlist((lapply(x, min)))), max(unlist((lapply(x, max)))), length=100),
                      contrasts = FALSE, nsim, ...) {
  if(!is.list(x) && length(x)<2) stop("At least two groups should be provided.")
  x.lengths <- as.numeric(lapply(x, FUN = length))
  if(!is.null(names(x))) groups <- rep(names(x), times=x.lengths)
  else groups <- rep(1:length(x), times=x.lengths)
  groups <- factor(groups, levels=unique(groups))
  gnames <- levels(groups)
  if(missing(nsim)) {
    if(!contrasts) {
      nsim <- length(x)*1000 - 1
    }
    else {
      J <- length(x)
      nsim <- (J*(J-1)/2)*1000 - 1
    }
    message("Creating ", nsim, " permutations.\n", sep="")
  }
  # setting the 'fun', "means" or "contrasts"
  if(!contrasts) fun <- ecdfmeans.m
  else fun <- ecdfcontrasts.m
  x <- unlist(x)
  # Observed difference between the ecdfs
  obs <- fun(x, groups, r)
  # Simulations by permuting to which groups each value belongs to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE), r), simplify = "array")
  complabels <- colnames(obs)

  csets <- list()
  for(i in 1:ncol(obs)) {
    csets[[i]] <- create_curve_set(list(r = r,
                                        obs = obs[,i],
                                        sim_m = sim[,i,]))
  }
  names(csets) <- complabels
  # GET
  res <- global_envelope_test(csets, alternative="two.sided", ..., nstep=1)
  if(!contrasts)
    res <- envelope_set_labs(res, xlab = expression(italic(x)),
                             ylab = expression(italic(hat(F)(x))))
  else
    res <- envelope_set_labs(res, xlab = expression(italic(x)),
                             ylab = expression(italic(hat(F)[i](x)-hat(F)[j](x))))
  attr(res, "contrasts") <- contrasts
  attr(res, "labels") <- complabels
  attr(res, "call") <- match.call()
  res
}
