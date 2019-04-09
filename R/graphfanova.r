# auxiliary functions
#--------------------

# statistics by group

# x: array(ndata x nt), each row is one functional data vector

vmean <- function(x) if (is.matrix(x)) apply(x, 2, mean) else x

vsum <- function(x) if (is.matrix(x)) apply(x, 2, sum) else x

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
  total <- apply(x^2, 2, sum) - apply(x, 2, sum)^2 / n
  within <- apply( groupSXX(x, groups) - groupSX(x, groups)^2/ni, 2, sum)
  (total / within - 1) / (k - 1) * (n-k)
}

#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats anova
corrFvalues <- function(x, groups) {
  Fvalues <- vector(length=ncol(x))
  for(i in 1:ncol(x)) {
    df <- data.frame(value = x[,i], group = groups)
    wl <- 1 / as.vector(by(df$value, df$group, function(x){ var(x, na.rm=T) }))
    df$w <- with(df, ifelse(group==1, wl[1], ifelse(group==2, wl[2], wl[3])))
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
#' This functions can be used to perform one-way graphical functional ANOVA tests described
#' in Mrkvička et al. (2016).
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
#' If you want to test the hypothesis
#' \deqn{H_0 : \mu_j(r) \equiv 0, j=1, \dots , J,}{H0: \mu_j(r) = 0, j=1,...,J,}
#' then you should use the test function
#' \deqn{\mathbf{T} = (\overline{T}_1({\bf r}), \overline{T}_2({\bf r}), \dots , \overline{T}_J({\bf r}))}{T = (\bar{T}_1(r), \bar{T}_2(r), ..., \bar{T}_J(r))}
#' where \eqn{\overline{T}_i({\bf r})}{\bar{T}_i(r)} is a vector of mean values of functions in the group j.
#' This can be done by choosing the summaryfun \code{"means"}.
#'
#' An alternative is to test the equivalent hypothesis
#' \deqn{H_0 : \mu_i(r) - \mu_j(r) = 0, i=1,\dots,J-1, j=1,\dots,J.}{H0: \mu_i(r) - \mu_j(r) = 0, i=1,...,J-1, j=i,...,J.}
#' This test corresponds to the post-hoc test done usually after an ANOVA test is significant, but
#' it can be directed tested by mean of the combined rank test (Mrkvička et al., 2017), if the
#' test vector is taken to consist of the differences of the group averages of test functions, namely
#' \deqn{\mathbf{T'} = (\overline{T}_1({\bf r})-\overline{T}_2({\bf r}),
#' \overline{T}_1({\bf r})-\overline{T}_3({\bf r}), \dots , \overline{T}_{J-1}({\bf r})-\overline{T}_J({\bf r})).}{T' = (\bar{T}_1(r)-\bar{T}_2(r), \bar{T}_1(r)-\bar{T}_3(r), ..., \bar{T}_{J-1}(r)-\bar{T}_J(r)).}
#' The summaryfun option \code{"contrasts"} can be used to perform the test based on this test vector.
#'
#' The test as such assumes that the variances are equal across the groups of functions. To deal with
#' unequal variances, the differences are rescaled as the first step as follows
#' \deqn{S_{ij}(r) = \frac{T_{ij}(r) - \overline{T}(r))}{\sqrt{\text{Var}(T_j(r))}} \sqrt{\text{Var}(T(r))} + \overline{T}(r))}{S_{ij}(r) = ( T_{ij}(r) - \bar{T}(r) ) / Sd(T_j(r)) * Sd(T(r)) + \bar{T}(r))}
#' where \eqn{\overline{T}({\bf r})}{\bar{T}(r)} is the overall sample mean and
#' \eqn{\sqrt{\text{Var}(T(r))}}{Sd(T(r))} is the overall sample standard deviation.
#' This scaling of the test functions can be obtained by giving the argument \code{variances = "unequal"}.
#'
#' @param nsim The number of random permutations.
#' @param curve_set The original data (an array of functions) provided as a \code{curve_set} object
#' (see \code{\link{create_curve_set}}) or a fdata object (see \code{\link[fda.usc]{fdata}}).
#' The curve set should include the argument values for the functions in the component \code{r}, and
#' the observed functions in the component \code{obs}.
#' @param groups The original groups (a factor vector representing the assignment to groups).
#' @param variances Either "equal" or "unequal". If "unequal", then correction for unequal variances
#' as explained in details will be done.
#' @param summaryfun Possible values are "means" and "contrasts".
#' See description for their meaning.
# Note: Possibly add a some arguments to specify which contrasts should be used.
# (Try to find our how this is usually done in R, in ordinary anova.)
#' @param n.aver If variances = "unequal", there is a possibility to use variances smoothed
#' by appying moving average to the estimated sample variances. n.aver determines
#' how many values on each side do contribute (incl. value itself).
#' @param mirror The complement of the argument circular of \code{\link[stats]{filter}}.
#' @param savefuns Logical. If TRUE, then the functions from permutations are saved to the attribute
#' simfuns.
#' @param test.equality A character with possible values \code{mean} (default), \code{var} and
#' \code{cov}. If \code{mean}, the functional ANOVA is performed to compare the means in the groups.
#' If \code{var}, then the equality of variances of the curves in the groups is tested by performing
#' the graphical functional ANOVA test on the functions
#' \deqn{Z_{ij}(r) = T_{ij}(r) - \bar{T}_j(r).}{Z_{ij}(r) = T_{ij}(r) - \bar{T}_j(r).}
#' If \code{cov}, then the equality of lag \code{cov.lag} covariance is tested by performing the fANOVA with
#' \deqn{W_{ij}(r) = \sqrt{|V_{ij}(r)|\cdot sign(V_{ij}(r))},}{|V_{ij}(r)| sign(V_{ij}(r)),}
#' where \deqn{V_{ij}(r) = (T_{ij}(r) - \bar{T}_j(r))((T_{ij}(r+s) - \bar{T}_j(r+s))).}{V_{ij}(r) = (T_{ij}(r) - \bar{T}_j(r))((T_{ij}(r+s) - \bar{T}_j(r+s))).}
#' See Mrkvicka et al. (2018) for more details.
#' @param cov.lag The lag of the covariance for testing the equality of covariances,
#' see \code{test.equality}.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @export
#' @seealso \code{\link{graph.fanova2d}}, \code{\link{frank.fanova}}
#' @references
#' Mrkvička, T., Hahn, U. and Myllymäki, M.
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#'
#' Mrkvička, T., Myllymäki, M., and Hahn, U. (2017).
#' Multiple Monte Carlo testing, with applications in spatial point processes.
#' Statistics and Computing 27 (5): 1239-1255. doi:10.1007/s11222-016-9683-9
#' @examples
#' #-- Centred government expenditure centralization ratios example
#' # This is an example analysis of the centred GEC in Mrkvicka et al.
#' data(cgec)
#' # Plot data in groups
#' subs <- function(group, ...) {
#'   cset <- cgec
#'   cset$obs <- cgec$obs[, attr(cset, "group") == group]
#'   plot(cset, ...)
#' }
#' par(mfrow=c(1,3))
#' for(i in 1:3) subs(i, main=paste("Group ", i, sep=""), ylab="Centred GEC")
#'
#' # Number of simulations
#' nsim <- 2499 # increase to reduce Monte Carlo error
#'
#' # Test for unequal lag 1 covariances
#' res.cov1 <- graph.fanova(nsim = nsim, curve_set = cgec,
#'                          groups = attr(cgec, "group"),
#'                          test.equality = "cov", cov.lag = 1)
#' plot(res.cov1, max_ncols_of_plots=3,
#'      labels = paste("Group ", 1:3, sep=""),
#'      xlab=substitute(paste(i, " (", italic(j), ")", sep=""), list(i="Year", j="r")),
#'      ylab=expression(italic(bar(W)[i](r))))
#' # Test for equality of variances among groups
#' res.var <- graph.fanova(nsim = nsim, curve_set = cgec,
#'                         groups = attr(cgec, "group"),
#'                         test.equality = "var")
#' plot(res.var, max_ncols_of_plots=3,
#'      labels = paste("Group ", 1:3, sep=""),
#'      xlab=substitute(paste(i, " (", italic(j), ")", sep=""), list(i="Year", j="r")),
#'      ylab=expression(italic(bar(Z)[i](r))))
#'
#' # Test for equality of means assuming equality of variances
#' # a) using 'means'
#' res <- graph.fanova(nsim = nsim, curve_set = cgec,
#'                     groups = attr(cgec, "group"),
#'                     variances = "equal",
#'                     summaryfun = "means")
#' plot(res, max_ncols_of_plots=3,
#'      labels = paste("Group ", 1:3, sep=""),
#'      xlab=substitute(paste(i, " (", italic(j), ")", sep=""), list(i="Year", j="r")),
#'      ylab=expression(italic(bar(T)[i](r))))
#' # b) using 'contrasts'
#' res2 <- graph.fanova(nsim = nsim, curve_set = cgec,
#'                     groups = attr(cgec, "group"),
#'                     variances = "equal",
#'                     summaryfun = "contrasts")
#' plot(res2, max_ncols_of_plots=3,
#'      xlab=substitute(paste(i, " (", italic(j), ")", sep=""), list(i="Year", j="r")),
#'      ylab=expression(italic(bar(T)[i](r)-bar(T)[j](r))))
#'
#' #-- Rimov water temperatures example
#' # This is an example analysis of the water temperature data set
#' # in Mrkvicka et al. (arXiv:1612.03608v2).
#' data(rimov)
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#' nsim <- 999
#'
#' # Test for equality of variances in the groups
#' resV <- graph.fanova(nsim=nsim, curve_set=rimov, groups=groups, summaryfun="means",
#'                      test.equality="var")
#' plot(resV)
#' # Test for equality of lag 1 covariances in the groups
#' resC <- graph.fanova(nsim=nsim, curve_set=rimov, groups=groups, summaryfun="means",
#'                      test.equality="cov", cov.lag=1)
#' plot(resC)
#'
#' # Test the equality of means in the groups (fANOVA), assuming equality of variances
#' res <- graph.fanova(nsim=nsim, curve_set=rimov, groups=groups, summaryfun="means")
#' plot(res)
#' res2 <- graph.fanova(nsim=nsim, curve_set=rimov, groups=groups, summaryfun="contrasts")
#' plot(res2)
graph.fanova <- function(nsim, curve_set, groups, variances="equal",
                         summaryfun = c("means", "contrasts"),
                         n.aver = 1L, mirror = FALSE, savefuns=FALSE,
                         test.equality = c("mean", "var", "cov"), cov.lag = 1, ...) {
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(!(class(curve_set) %in% c("curve_set", "fdata"))) stop("The curve_set does not have a valid class.\n")
  curve_set <- convert_fdata(curve_set)
  if(!is.matrix(curve_set[['obs']])) stop("The curve_set must include data functions (sim_m ignored).\n")
  x <- t(curve_set[['obs']])
  if(nrow(x) != length(groups)) stop("The length of groups should be equal with the number of functions.\n")
  if(!is.factor(groups)) {
    warning("The argument groups is not a factor. Transforming it to a factor by as.factor.\n")
    groups <- as.factor(groups)
  }
  r <- curve_set[['r']]
  test.equality <- match.arg(test.equality)
  switch(test.equality,
         "mean" = {
           if(!(variances %in% c("equal", "unequal"))) stop("Options for variances are equal and unequal.\n")
           if(variances == "unequal") x <- corrUnequalVar(x, groups, n.aver, mirror)
         },
         "var" = {
           x <- testUnequalVarTrans(x, groups)
         },
         "cov" = {
           x <- testUnequalCovTrans(x, groups, lag=cov.lag)
           r <- r[1:(length(r)-cov.lag)]
         })

  summaryfun <- match.arg(summaryfun)
  # setting that 'summaryfun' is a function
  switch(summaryfun, 
         means = {fun = means.m},
         contrasts = {fun = contrasts.m}
         )

  obs <- fun(x, groups)
  # simulations by permuting to which groups the functions belong to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE)), simplify = "array")
  # labels for comparisons (rownames(sim) is the same as names(obs))
  complabels <- colnames(obs)

  csets <- NULL
  for(i in 1:ncol(obs)) {
    csets[[i]] <- create_curve_set(list(r = r,
                                        obs = obs[,i],
                                        sim_m = sim[,i,]))
  }
  res <- global_envelope_test(csets, alternative="two.sided", ..., nstep=1)
  if(inherits(res, "combined_global_envelope")) names(res) <- complabels
  attr(res, "method") <- "Graphical functional ANOVA" # Change method name
  attr(res, "summaryfun") <- summaryfun
  attr(res, "labels") <- complabels
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
#' \eqn{n_1,\dosts,n_J}{n1, ..., nJ} functions
#' \eqn{T_{ij}, i=\dots,J, j=1,\dots,n_j}{T_{ij}, i=1,...,J, j=1,...,nj}.
#' The functions should be given in the argument x, and the groups in the argument groups.
#' The test assumes that there exists non random functions \eqn{\mu(r)}{\mu(r)} and
#' \eqn{\mu_i(r)}{\mu_i(r)} such that
#' \deqn{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, \dots, J, j=1, \dots , n_j}{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, ..., J, j=1, ..., nj}
#' where \eqn{e_{ij}(r)}{e_{ij}(r)} are independent and normally distributed.
#' The test vector is
#' \deqn{\mathbf{T} = (F(r_1), F(r_2), \dots , F(r_K))}{T = (F(r_1), F(r_2), \dots , F(r_K))},
#' where \eqn{F(r_i)}{F(r_i)} stands for the F-statistic. The simulations are performed by
#' permuting the test functions. Further details can be found in Mrkvička et al. (2016).
#'
#' The argument \code{equalvar=TRUE} means that equal variances across groups are assumed.
#' The correction for unequal variances can be done by using the corrected F-statistic
#' (option \code{equalvar=FALSE}).
#'
#' Unfortunately this test is not able to detect which groups are different from each other.
#'
#' @inheritParams graph.fanova
#' @param variances Either "equal" or "unequal". If "equal", then the traditional F-values are used.
#' If "unequal", then the corrected F-values are used. The current implementation uses
#' \code{\link[stats]{lm}} to get the corrected F-values.
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2018)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#' @examples
#' \donttest{
#' data(rimov)
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#' res <- frank.fanova(nsim=2499, curve_set=rimov, groups=groups)
#' plot(res, ylab="F-statistic")
#' }
frank.fanova <- function(nsim, curve_set, groups, variances="equal",
                         test.equality = c("mean", "var", "cov"), cov.lag = 1, ...) {
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(!(class(curve_set) %in% c("curve_set", "fdata"))) stop("The curve_set does not have a valid class.\n")
  curve_set <- convert_fdata(curve_set)
  if(!is.matrix(curve_set[['obs']])) stop("The curve_set must include data functions (sim_m ignored).\n")
  x <- t(curve_set[['obs']])
  if(nrow(x) != length(groups)) stop("The length of groups should be equal with the number of functions.\n")
  if(!is.factor(groups)) {
    warning("The argument groups is not a factor. Transforming it to a factor by as.factor.\n")
    groups <- as.factor(groups)
  }
  r <- curve_set[['r']]
  test.equality <- match.arg(test.equality)
  switch(test.equality,
         "mean" = {
           if(!(variances %in% c("equal", "unequal"))) stop("Options for variances are equal and unequal.\n")
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
  global_envelope_test(cset, alternative="greater", ...)
}

#' One-way graphical functional ANOVA for images
#'
#' One-way ANOVA tests for image data with graphical interpretation
#'
#'
#' This function can be used to perform one-way graphical functional ANOVA tests described
#' in Mrkvička et al. (2018). The function transforms images to vectors (1d), calls
#' \code{\link{graph.fanova}} and transform results back to images (2d).
#'
#' @inheritParams graph.fanova
#' @inheritParams frank.fanova2d
#' @param ... Additional parameters to be passed to \code{\link{graph.fanova}}.
#' @return A \code{global_envelope_2d} or \code{combined_global_envelope_2d} object, which can be
#' printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2018)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#' @examples
#' data("imageset1")
#' res <- graph.fanova2d(nsim = 19, # Increase nsim for serious analysis!
#'                       image_set = imageset1$image_set,
#'                       groups = imageset1$Group)
#' par(mfrow=c(2,5))
#' plot(res)
#' # Contrasts
#' res.c <- graph.fanova2d(nsim = 19, # Increase nsim for serious analysis!
#'                         image_set = imageset1$image_set,
#'                         groups = imageset1$Group,
#'                         summaryfun = "contrasts")
#' par(mfrow=c(2,3))
#' plot(res.c)
graph.fanova2d <- function(nsim, image_set, groups, ...) {
  if(class(image_set)[1] != "image_set") stop("The image_set does not have a valid class.\n")
  obs_d <- dim(image_set$obs)
  sim_d <- dim(image_set$sim_m)
  image_set <- check_image_set_dimensions(image_set)
  # Create curve sets transforming the 2d functions (matrices) to vectors
  curve_set_v <- image_set_to_curve_set(image_set)

  args <- list(...)
  if("test.equality" %in% names(args))
    if(!(args$test.equality %in% c("mean", "var"))) stop("Only options \'mean\' and \'var\' available for test.equality for images.\n")
  res_v <- graph.fanova(nsim=nsim, curve_set=curve_set_v, groups=groups, ..., n.aver=1)
  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_set)
  attr(res, "call") <- match.call()
  res
}

#' Rank envelope F-test for images
#'
#' Rank envelope F-test for images
#'
#'
#' A one-way functional ANOVA for images (2d functions), based on the rank envelope applied to
#' F values.
#' This function can be used to perform F-rank one-way functional ANOVA tests described
#' in Mrkvička et al. (2018). The function transforms images to vectors (1d), calls
#' \code{\link{frank.fanova}} and transform results back to images (2d).
#'
#' @inheritParams frank.fanova
#' @param image_set A set of images containing the data, see \code{\link{create_image_set}}.
#' @param ... Additional parameters to be passed to \code{\link{frank.fanova}}.
#' @export
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2018)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#' @examples
#' data("imageset1")
#' res <- frank.fanova2d(nsim = 19, # Increase nsim for serious analysis!
#'                       image_set = imageset1$image_set,
#'                       groups = imageset1$Group)
#' par(mfrow=c(1,3))
#' plot(res)
frank.fanova2d <- function(nsim, image_set, groups, ...) {
  if(class(image_set)[1] != "image_set") stop("The image_set does not have a valid class.\n")
  obs_d <- dim(image_set$obs)
  image_set <- check_image_set_dimensions(image_set)
  # Create curve sets transforming the 2d functions (matrices) to vectors
  curve_set_v <- image_set_to_curve_set(image_set)

  res_v <- frank.fanova(nsim=nsim, curve_set=curve_set_v, groups=groups, ...)
  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_set)
  attr(res, "call") <- match.call()
  res
}
