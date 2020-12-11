# auxiliary functions
#--------------------

# statistics by group

# x: array(ndata x nt), each row is one functional data vector

vmean <- function(x) if (is.matrix(x)) colMeans(x) else x

vsum <- function(x) if (is.matrix(x)) colSums(x) else x

#' @importFrom stats var
vvar <- function(x) if (is.matrix(x)) apply(x, 2, stats::var) else 0*x

groupmeans <- function(x, groups) t(sapply(levels(groups), function(g) vmean(x[groups==g,])))

groupvar <- function(x, groups) t(sapply(levels(groups), function(g) vvar(x[groups==g,])))

groupn <- function(groups) sapply(levels(groups), function(g) sum(groups == g))

groupSX <- function(x, groups) t(sapply(levels(groups), function(g) vsum(x[groups==g,])))

groupSXX <- function(x, groups) t(sapply(levels(groups), function(g) vsum(x[groups==g,]^2)))

grouperror <- function(x, groups) groupvar(x, groups) / groupn(groups)

# Moving average, for variance estimation.
#
# @param x
# @param naver How many on each side do contribute (incl value itself)
# @param mirror Mirror at boundaries. If TRUE, then circular, i.e. wrap the filter around the ends of the functions.
#' @importFrom stats filter
maverage <- function(x, n.aver = 1L, mirror = FALSE) {
  if(n.aver >= (length(x)+1)%/%2) return(rep(mean(x), length(x)))
  if(n.aver == 1) return(x)
  lfilt <- n.aver * 2 - 1
  if(mirror) x <- c(x[n.aver:2], x, x[(-1:-(n.aver -1))+length(x)])
  mav <- stats::filter(x, rep(1, lfilt)/lfilt, method = "convolution", sides = 2, circular = !mirror)
  if(mirror) mav <- mav[-c(1:(n.aver-1), length(mav)+(0: (2-n.aver)))]
  mav
}

# Transformation to equalize variances in groups
corrUnequalVar <- function(x, groups, ...) {
  # Group means
  m <- groupmeans(x, groups)
  # Sample variance over all functions, Var(T(r))
  varT <- vvar(x)
  # Variances in the groups, Var(T_j(r))
  v <- groupvar(x, groups)
  # Moving average
  averargs <- list(...)
  if(length(averargs) > 0){
    v <- t(apply(v, 1, maverage, ...))
    varT <- maverage(varT, ...)
  }
  # Take S_ij(r) = (T_ij(r) - \bar{T}_j(r)) / \sqrt( Var(T_j(r)) ) * Var(T(r)) + \bar{T}_j(r)
  for(i in 1:nrow(x)) {
    x[i,] <- (x[i,] - m[which(rownames(m) == groups[i]),]) / sqrt(v[which(rownames(v) == groups[i]),]) * sqrt(varT) + m[which(rownames(m) == groups[i]),]
  }
  x
}

# Transformation for testing equality of variances in groups
# Take Z_ij(r) = |T_ij(r) - \bar{T}_j(r))|
testUnequalVarTrans <- function(x, groups) {
  # Group means
  m <- groupmeans(x, groups)
  # Take Z_ij(r) = |T_ij(r) - \bar{T}_j(r))|
  for(i in 1:nrow(x)) {
    x[i,] <- abs(x[i,] - m[which(rownames(m) == groups[i]),])
  }
  x
}

# Tranformation for testing equality of lag s covariances in groups
testUnequalCovTrans <- function(x, groups, lag=1) {
  n.obs <- ncol(x)
  # Group means
  m <- groupmeans(x, groups)
  # Take Z_ij(r) = sqrt{(T_ij(r) - \bar{T}_j(r)))(T_ij(r+lag) - \bar{T}_j(r+lag)))} * sign
  xnew <- matrix(nrow=nrow(x), ncol=ncol(x)-lag)
  for(i in 1:nrow(xnew)) {
    xj <- m[which(rownames(m) == groups[i]),] # Group mean
    cov <- (x[i,1:(n.obs-lag)] - xj[1:(n.obs-lag)])*(x[i,(1+lag):n.obs]-xj[(1+lag):n.obs])
    xnew[i,] <- sqrt(abs(cov)) * sign(cov)
  }
  xnew
}

# group statistics
#-----------------
# x = An array with the original functions
# groups = a factor vector representing the assignment to groups
# ... Ignored

Fvalues <- function(x, groups) {
  ni <- groupn(groups)
  n <- sum(ni)
  k <- nlevels(groups)
  total <- colSums(x^2) - colSums(x)^2 / n
  within <- colSums( groupSXX(x, groups) - groupSX(x, groups)^2/ni )
  (total / within - 1) / (k - 1) * (n-k)
}

#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats anova
corrFvalues <- function(x, groups) {
  Fvalues <- vector(length=ncol(x))
  lgroups <- levels(groups)
  for(i in 1:ncol(x)) {
    df <- data.frame(value = x[,i], group = groups)
    wl <- 1 / as.vector(by(df$value, df$group, function(x){ var(x, na.rm=T) }))
    df$w <- wl[sapply(1:nrow(df), FUN = function(i) { which(lgroups %in% df$group[i]) })]
    w.mod <- stats::lm(value~group, data=df, na.action=stats::na.omit, weights=df$w)
    temp <- stats::anova(w.mod)
    Fvalues[i] <- temp$F[1] # anova(aov(formula = Lvalues ~ group, data = df))$F[1]
  }
  Fvalues
}

#- means and contrasts as matrices

means.m <- function(x, groups, ...) {
  jm <- t(groupmeans(x, groups))
  colnames(jm) <- levels(groups)
  jm
}

contrasts.m <- function(x, groups, ...) {
  k <- nlevels(groups)
  gnam <- levels(groups)
  mea <- groupmeans(x, groups)
  cont <- matrix(0, nrow=dim(x)[2], ncol=k*(k-1)/2)
  cont.names <- vector(length=k*(k-1)/2)
  counter <- 1
  for(i in 1:(k-1)) for(j in (i+1):k) {
    cont[, counter] <- mea[i, ] - mea[j, ]
    cont.names[counter] <- paste(gnam[i], gnam[j], sep="-")
    counter <- counter+1
  }
  colnames(cont) <- cont.names
  cont
}


#' One-way graphical functional ANOVA
#'
#' One-way ANOVA tests for functional data with graphical interpretation
#'
#'
#' This function can be used to perform one-way graphical functional ANOVA tests described
#' in Mrkvička et al. (2020). Both 1d and 2d functions are allowed in curve sets.
#'
#' The tests assume that there are \eqn{J}{J} groups which contain
#' \eqn{n_1,\dots,n_J}{n1, ..., nJ} functions
#' \eqn{T_{ij}, i=\dots,J, j=1,\dots,n_j}{T_{ij}, i=1,...,J, j=1,...,nj}.
#' The functions should be given in the argument \code{curve_set},
#' and the groups in the argument \code{groups}.
#' The tests assume that \eqn{T_{ij}, i=1,...,n_j}{T_{ij}, i=1,...,n_j} is an iid sample from
#' a stochastic process with mean function \eqn{\mu_j}{\mu_j} and
#' covariance function \eqn{\gamma_j(s,t)}{\gamma_j(s,t)} for s,t in R and j = 1,..., J.
#'
#' To test the hypothesis
#' \deqn{H_0 : \mu_1(r) \equiv \mu_2(r)\equiv \dots \equiv \mu_J(r),}{H0: \mu_1(r) = \mu_2(r) = ... = \mu_J(r),}
#' you can use the test function
#' \deqn{\mathbf{T} = (\overline{T}_1({\bf r}), \overline{T}_2({\bf r}), \dots , \overline{T}_J({\bf r}))}{T = (\bar{T}_1(r), \bar{T}_2(r), ..., \bar{T}_J(r))}
#' where \eqn{\overline{T}_i({\bf r})}{\bar{T}_i(r)} is a vector of mean values of functions in the group j.
#' This test function is used when \code{contrasts = FALSE} (default).
#'
#' The hypothesis can equivalently be written as
#' \deqn{H_0 : \mu_i(r) - \mu_j(r) = 0, i=1,\dots,J-1, j=1,\dots,J.}{H0: \mu_i(r) - \mu_j(r) = 0, i=1,...,J-1, j=i,...,J.}
#' and, alternatively, one can use the test function (vector)
#' taken to consist of the differences of the group averages,
#' \deqn{\mathbf{T'} = (\overline{T}_1({\bf r})-\overline{T}_2({\bf r}),
#' \overline{T}_1({\bf r})-\overline{T}_3({\bf r}), \dots , \overline{T}_{J-1}({\bf r})-\overline{T}_J({\bf r})).}{T' = (\bar{T}_1(r)-\bar{T}_2(r), \bar{T}_1(r)-\bar{T}_3(r), ..., \bar{T}_{J-1}(r)-\bar{T}_J(r)).}
#' The choice is available with the option \code{contrasts = TRUE}.
#' This test corresponds to the post-hoc test done usually after an ANOVA test is significant, but
#' it can be directed tested by means of the combined rank test (Mrkvička et al., 2017) with this test vector.
#'
#' The test as such assumes that the variances are equal across the groups of functions. To deal with
#' unequal variances, the differences are rescaled as the first step as follows
#' \deqn{S_{ij}(r) = \frac{T_{ij}(r) - \overline{T}(r))}{\sqrt{Var(T_j(r))}} \sqrt{Var(T(r))} + \overline{T}(r))}{S_{ij}(r) = ( T_{ij}(r) - \bar{T}(r) ) / Sd(T_j(r)) * Sd(T(r)) + \bar{T}(r))}
#' where \eqn{\overline{T}({\bf r})}{\bar{T}(r)} is the overall sample mean and
#' \eqn{\sqrt{Var(T(r))}}{Sd(T(r))} is the overall sample standard deviation.
#' This scaling of the test functions can be obtained by giving the argument \code{variances = "unequal"}.
#'
#' @param nsim The number of random permutations.
#' @param curve_set The original data (an array of functions) provided as a \code{curve_set} object
#' (see \code{\link{create_curve_set}}) or a fdata object (see \code{\link[fda.usc]{fdata}}).
#' The curve set should include the argument values for the functions in the component \code{r}, and
#' the observed functions in the component \code{obs}.
#' @param groups The original groups (a factor vector representing the assignment to groups).
#' @param variances Either "equal" or "unequal". If "unequal", then correction for unequal variances
#' as explained in details will be done. Only relevant for the case \code{test.equality = "means"} (default).
#' @param contrasts Logical. FALSE and TRUE specify the two test functions as described in
#' description part of this help file.
# Note: Possibly add a some arguments to specify which contrasts should be used.
# (Try to find our how this is usually done in R, in ordinary anova.)
#' @param n.aver If variances = "unequal", there is a possibility to use variances smoothed
#' by appying moving average to the estimated sample variances. n.aver determines
#' how many values on each side do contribute (incl. value itself).
#' @param mirror The complement of the argument circular of \code{\link[stats]{filter}}. Another parameter
#' for the moving average to estimate sample variances (see \code{n.aver}).
#' @param savefuns Logical. If TRUE, then the functions from permutations are saved to the attribute
#' simfuns.
#' @param test.equality A character with possible values \code{mean} (default), \code{var} and
#' \code{cov}. If \code{mean}, the functional ANOVA is performed to compare the means in the groups.
#' If \code{var}, then the equality of variances of the curves in the groups is tested by performing
#' the graphical functional ANOVA test on the functions
#' \deqn{Z_{ij}(r) = T_{ij}(r) - \bar{T}_j(r).}{Z_{ij}(r) = T_{ij}(r) - \bar{T}_j(r).}
#' If \code{cov}, then the equality of lag \code{cov.lag} covariance is tested by performing the fANOVA with
#' \deqn{W_{ij}(r) = \sqrt{|V_{ij}(r)|\cdot sign(V_{ij}(r))},}{|V_{ij}(r)| sign(V_{ij}(r)),}
#' where
#' \deqn{V_{ij}(r) = (T_{ij}(r) - \bar{T}_j(r))((T_{ij}(r+s) - \bar{T}_j(r+s))).}{V_{ij}(r) = (T_{ij}(r) - \bar{T}_j(r))((T_{ij}(r+s) - \bar{T}_j(r+s))).}
#' See Mrkvicka et al. (2020) for more details.
#' @param cov.lag The lag of the covariance for testing the equality of covariances,
#' see \code{test.equality}.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @export
#' @seealso \code{\link{frank.fanova}}
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#'
#' Mrkvička, T., Myllymäki, M., and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics and Computing 27 (5): 1239-1255. doi:10.1007/s11222-016-9683-9
#'
#' Myllymäki, M and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#' @examples
#' #-- NOx levels example (see for details Myllymaki and Mrkvicka, 2020)
#' if(require("fda.usc", quietly=TRUE)) {
#'   # Prepare data
#'   data("poblenou")
#'   Free <- poblenou$df$day.festive == 1 |
#'     as.integer(poblenou$df$day.week) >= 6
#'   MonThu <- poblenou$df$day.festive == 0 & poblenou$df$day.week %in% 1:4
#'   Friday <- poblenou$df$day.festive == 0 & poblenou$df$day.week == 5
#'   Type <- vector(length=length(Free))
#'   Type[Free] <- "Free"
#'   Type[MonThu] <- "MonThu"
#'   Type[Friday] <- "Fri"
#'   Type <- factor(Type, levels = c("MonThu", "Fri", "Free"))
#'
#'   # (log) Data as a curve_set
#'   cset <- create_curve_set(list(r=0:23,
#'              obs=t(log(poblenou[['nox']][['data']]))))
#'   # Graphical functional ANOVA
#'  \dontshow{nsim <- 19}
#'  \donttest{nsim <- 2999}
#'   res.c <- graph.fanova(nsim = nsim, curve_set = cset,
#'                         groups = Type, variances = "unequal",
#'                         contrasts = TRUE)
#'   plot(res.c) + ggplot2::labs(x = "Hour", y = "Diff.")
#' }
#'
#' #-- Centred government expenditure centralization ratios example
#' # This is an example analysis of the centred GEC in Mrkvicka et al.
#' data("cgec")
#'
#' # Number of simulations
#' \dontshow{nsim <- 19}
#' \donttest{nsim <- 2499 # increase to reduce Monte Carlo error}
#'
#' # Test for unequal lag 1 covariances
#' res.cov1 <- graph.fanova(nsim = nsim, curve_set = cgec$cgec,
#'                          groups = cgec$group,
#'                          test.equality = "cov", cov.lag = 1)
#' plot(res.cov1)
#' # Add labels
#' plot(res.cov1, labels = paste("Group ", 1:3, sep="")) +
#'   ggplot2::xlab(substitute(paste(italic(i), " (", j, ")", sep=""), list(i="r", j="Year")))
#' # Test for equality of variances among groups
#' res.var <- graph.fanova(nsim = nsim, curve_set = cgec$cgec,
#'                         groups = cgec$group,
#'                         test.equality = "var")
#' plot(res.var)
#'
#' # Test for equality of means assuming equality of variances
#' # a) using 'means'
#' res <- graph.fanova(nsim = nsim, curve_set = cgec$cgec,
#'                     groups = cgec$group,
#'                     variances = "equal",
#'                     contrasts = FALSE)
#' plot(res)
#' # b) using 'contrasts'
#' res2 <- graph.fanova(nsim = nsim, curve_set = cgec$cgec,
#'                      groups = cgec$group,
#'                      variances = "equal",
#'                      contrasts = TRUE)
#' plot(res2)
#'
#' # Image set examples
#' data("imageset3")
#' res <- graph.fanova(nsim = 19, # Increase nsim for serious analysis!
#'                     curve_set = imageset3$image_set,
#'                     groups = imageset3$Group)
#' plot(res)
#' # Contrasts
#' res.c <- graph.fanova(nsim = 19, # Increase nsim for serious analysis!
#'                       curve_set = imageset3$image_set,
#'                       groups = imageset3$Group,
#'                       contrasts = TRUE)
#' plot(res.c)
graph.fanova <- function(nsim, curve_set, groups, variances="equal",
                         contrasts = FALSE,
                         n.aver = 1L, mirror = FALSE, savefuns=FALSE,
                         test.equality = c("mean", "var", "cov"), cov.lag = 1, ...) {
  if(nsim < 1) stop("Not a reasonable value of nsim.")
  if(!inherits(curve_set, c("curve_set", "fdata"))) stop("The curve_set does not have a valid class.")
  curve_set <- convert_fdata(curve_set)
  if(curve_set_is1obs(curve_set)) stop("All (data) functions of the curve_set must be equal.")
  x <- data_and_sim_curves(curve_set)
  if(nrow(x) != length(groups)) stop("The length of groups should be equal with the number of functions.")
  if(!is.factor(groups)) {
    warning("The argument groups is not a factor. Transforming it to a factor by as.factor.")
    groups <- as.factor(groups)
  }
  test.equality <- match.arg(test.equality)
  if(test.equality=="cov" && !curve_set_is1d(curve_set)) {
    stop("test.equality='cov' is only available for 1 dimensional functions")
  }
  r <- curve_set[['r']]
  switch(test.equality,
         "mean" = {
           if(!(variances %in% c("equal", "unequal"))) stop("Options for variances are equal and unequal.")
           if(variances == "unequal") x <- corrUnequalVar(x, groups, n.aver, mirror)
           if(!contrasts) ylab <- expression(italic(bar(T)[i](r)))
           else ylab <- expression(italic(bar(T)[i](r)-bar(T)[j](r)))
         },
         "var" = {
           x <- testUnequalVarTrans(x, groups)
           if(!contrasts) ylab <- expression(italic(bar(Z)[i](r)))
           else ylab <- expression(italic(bar(Z)[i](r)-bar(Z)[j](r)))
         },
         "cov" = {
           x <- testUnequalCovTrans(x, groups, lag=cov.lag)
           r <- r[1:(length(r)-cov.lag)]
           if(!contrasts) ylab <- expression(italic(bar(W)[i](r)))
           else ylab <- expression(italic(bar(W)[i](r)-bar(W)[j](r)))
         })

  # setting the 'fun', "means" or "constrasts"
  if(!contrasts) fun <- means.m
  else fun <- contrasts.m

  obs <- fun(x, groups)
  # simulations by permuting to which groups the functions belong to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE)), simplify = "array")
  # labels for comparisons (rownames(sim) is the same as names(obs))
  complabels <- colnames(obs)

  csets <- NULL
  for(i in 1:ncol(obs)) {
    csets[[complabels[i]]] <- create_curve_set(list(r = r,
                                        obs = obs[,i],
                                        sim_m = sim[,i,]))
  }
  res <- global_envelope_test(csets, alternative="two.sided", ..., nstep=1)
  attr(res, "method") <- "Graphical functional ANOVA" # Change method name
  attr(res, "contrasts") <- contrasts
  attr(res, "labels") <- complabels
  res <- envelope_set_labs(res, ylab=ylab)
  attr(res, "call") <- match.call()
  if(savefuns) attr(res, "simfuns") <- csets
  res
}

#' Rank envelope F-test
#'
#' A one-way functional ANOVA based on the rank envelope applied to F values
#'
#'
#' The test assumes that there are \eqn{J}{J} groups which contain
#' \eqn{n_1,\dots,n_J}{n1, ..., nJ} functions
#' \eqn{T_{ij}, i=\dots,J, j=1,\dots,n_j}{T_{ij}, i=1,...,J, j=1,...,nj}.
#' The functions should be given in the argument x, and the groups in the argument groups.
#' The test assumes that there exists non random functions \eqn{\mu(r)}{\mu(r)} and
#' \eqn{\mu_i(r)}{\mu_i(r)} such that
#' \deqn{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, \dots, J, j=1, \dots , n_j}{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, ..., J, j=1, ..., nj}
#' where \eqn{e_{ij}(r)}{e_{ij}(r)} are independent and normally distributed.
#' The test vector is
#' \deqn{\mathbf{T} = (F(r_1), F(r_2), \dots , F(r_K)),}{T = (F(r_1), F(r_2), \dots , F(r_K)),}
#' where \eqn{F(r_i)}{F(r_i)} stands for the F-statistic. The simulations are performed by
#' permuting the test functions. Further details can be found in Mrkvička et al. (2020).
#'
#' The argument \code{variances="equal"} means that equal variances across groups are assumed.
#' The correction for unequal variances can be done by using the corrected F-statistic
#' (option \code{variances="unequal"}).
#'
#' Unfortunately this test is not able to detect which groups are different from each other.
#'
#' @inheritParams graph.fanova
#' @param variances Either "equal" or "unequal". If "equal", then the traditional F-values are used.
#' If "unequal", then the corrected F-values are used. The current implementation uses
#' \code{\link[stats]{lm}} to get the corrected F-values.
#' @seealso graph.fanova
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#' @examples
#' data("rimov")
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#' \donttest{res <- frank.fanova(nsim=2499, curve_set=rimov, groups=groups)}
#' \dontshow{res <- frank.fanova(nsim=4, curve_set=rimov, groups=groups, alpha=0.2)}
#' plot(res, ylab="F-statistic")
#'
#' data("imageset3")
#' res2 <- frank.fanova(nsim = 19, # Increase nsim for serious analysis!
#'                      curve_set = imageset3$image_set,
#'                      groups = imageset3$Group)
#' plot(res2)
#' plot(res2, fixedscales=FALSE)
frank.fanova <- function(nsim, curve_set, groups, variances="equal",
                         test.equality = c("mean", "var", "cov"), cov.lag = 1, ...) {
  if(nsim < 1) stop("Not a reasonable value of nsim.")
  if(!inherits(curve_set, c("curve_set", "fdata"))) stop("The curve_set does not have a valid class.")
  curve_set <- convert_fdata(curve_set)
  if(curve_set_is1obs(curve_set)) stop("All (data) functions of the curve_set must be equal.")
  x <- data_and_sim_curves(curve_set)
  if(nrow(x) != length(groups)) stop("The length of groups should be equal with the number of functions.")
  if(!is.factor(groups)) {
    warning("The argument groups is not a factor. Transforming it to a factor by as.factor.")
    groups <- as.factor(groups)
  }
  test.equality <- match.arg(test.equality)
  if(test.equality=="cov" && !curve_set_is1d(curve_set)) {
    stop("test.equality='cov' is only available for 1 dimensional functions")
  }
  r <- curve_set[['r']]
  switch(test.equality,
         "mean" = {
           if(!(variances %in% c("equal", "unequal"))) stop("Options for variances are equal and unequal.")
           if(variances == "equal") fun <- Fvalues
           else fun <- corrFvalues
         },
         "var" = {
           x <- testUnequalVarTrans(x, groups)
           fun <- Fvalues
         },
         "cov" = {
           x <- testUnequalCovTrans(x, groups, lag=cov.lag)
           r <- r[1:(length(r)-cov.lag)]
           fun <- Fvalues
         })

  obs <- fun(x, groups)
  # simulations by permuting to which groups the functions belong to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE)))

  cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))
  # Perform the global envelope test
  res <- global_envelope_test(cset, alternative="greater", ...)
  res <- envelope_set_labs(res, ylab = expression(italic(F(r))))
  res
}
