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

# group statistics
#-----------------
# x = An array with the original functions
# groups = a factor vector representing the assignment to groups

#' @importFrom fda.usc is.fdata
Fvalues <- function(x, groups) {
  if (fda.usc::is.fdata(x)) x <- x$data
  ni <- groupn(groups)
  n <- sum(ni)
  k <- nlevels(groups)
  total <- apply(x^2, 2, sum) - apply(x, 2, sum)^2 / n
  within <- apply( groupSXX(x, groups) - groupSX(x, groups)^2/ni, 2, sum)
  (total / within - 1) / (k - 1) * (n-k)
}

#' @importFrom fda.usc is.fdata
#' @importFrom stats lm
corrFvalues <- function(x, groups) {
  if (fda.usc::is.fdata(x)) x <- x$data
  Fvalues <- vector(length=ncol(x))
  for(i in 1:ncol(x)) {
    df <- data.frame(value = x[,i], group = groups)
    wl <- 1 / as.vector(by(df$value, df$group, function(x){ var(x, na.rm=T) }))
    df$w <- with(df, ifelse(group==1, wl[1], ifelse(group==2, wl[2], wl[3])))
    w.mod <- stats::lm(value~group, df, na.action=na.omit, weights=w)
    temp <- anova(w.mod)
    Fvalues[i] <- temp$F[1] # anova(aov(formula = Lvalues ~ group, data = df))$F[1]
  }
  Fvalues
}


#- means as long vector

# ... Ignored
#' @importFrom fda.usc is.fdata
means <- function(x, groups, ...){
  if(fda.usc::is.fdata(x)) x <- x$data
  jm <- as.vector(t(groupmeans(x, groups)))
  names(jm) <- rep(levels(groups), each = dim(x)[2])
  jm
}

#' @importFrom fda.usc is.fdata
studmeans <- function(x, groups, ...){
  if(fda.usc::is.fdata(x)) x <- x$data
  # \bar{T}_i
  jm <- as.vector(t(groupmeans(x, groups)))
  names(jm) <- rep(levels(groups), each = dim(x)[2])
  # Var(\bar{T}_i)
  err <- grouperror(x, groups)
  glevels <- levels(groups)
  # Var(\bar{T}_{-i}(r))
  err.minusi <- t(sapply(levels(groups), function(g) vvar(x[groups!=g,]) / sum(groups != g) ))
  # Moving average
  averargs <- list(...)
  if(length(averargs) > 0){
    err <- t(apply(err, 1, maverage, ...))
    err.minusi <- t(apply(err.minusi, 1, maverage, ...))
  }
  # ( \bar{T}_i - \bar{T}_{-i}(r) ) / (Var(\bar{T}_i) + Var(\bar{T}_{-i}(r)))
  for(i in 1:length(glevels)) {
    T.minusi <- apply(x[groups != glevels[i],], MARGIN=2, FUN=mean) # mean \bar{T}_{-i}(r)
    jm[names(jm) == glevels[i]] <- ( jm[names(jm) == glevels[i]] - T.minusi ) / sqrt(err[i,] + err.minusi[i,])
  }
  jm
}


#- contrasts as long vectors

# ... Ignored
#' @importFrom fda.usc is.fdata
contrasts <- function(x, groups, ...){
  if (fda.usc::is.fdata(x)) x <- x$data
  k <- nlevels(groups)
  gnam <- levels(groups)
  mea <- groupmeans(x, groups)
  cont <- NULL
  nt <- dim(x)[2]
  for(i in 1:(k-1)) for(j in (i+1):k)  {
    ct <- mea[i, ] - mea[j, ]
    names(ct) <- rep(paste(gnam[i], gnam[j], sep="-"), nt)
    cont <- c(cont, ct)
  }
  cont
}

# @param ... Additional parameters passed to maverage
#' @importFrom fda.usc is.fdata
#' @importFrom caTools runmean
studcontrasts <- function(x, groups, ...){
  if (fda.usc::is.fdata(x)) x <- x$data
  k <- nlevels(groups)
  gnam <- levels(groups)
  mea <- groupmeans(x, groups)
  err <- grouperror(x, groups)
  averargs <- list(...)
  if(length(averargs) > 0){
    err <- t(apply(err, 1, maverage, ...))
  }
  nt <- dim(x)[2]
  scont <- NULL
  for(i in 1:(k-1)) for (j in (i+1):k){
    sct <- (mea[i, ] - mea[j, ])/ sqrt(err[i]+err[j])
    names(sct) <- rep(paste(gnam[i], gnam[j], sep="-"), nt)
    scont <- c(scont, sct)
  }
  scont
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
#' \eqn{n_1,\dosts,n_J}{n1, ..., nJ} functions
#' \eqn{T_{ij}, i=\dots,J, j=1,\dots,n_j}{T_{ij}, i=1,...,J, j=1,...,nj}.
#' The functions should be given in the argument x, and the groups in the argument groups.
#' The tests assume that there exists non random functions \eqn{\mu(r)}{\mu(r)} and
#' \eqn{\mu_i(r)}{\mu_i(r)} such that
#' \deqn{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, \dots, J, j=1, \dots , n_j}{T_{ij}(r) =\mu(r) + \mu_i(r) + e_{ij}(r), i=1, ..., J, j=1, ..., nj}
#' where \eqn{e_{ij}(r)}{e_{ij}(r)} are random errors (either correlated or not).
#'
#' If you want to test the hypothesis
#' \deqn{H_0 : \mu_i(r) \equiv 0, i=1, \dots , J,}{H0: \mu_i(r) = 0, i=1,...,J,}
#' then you should use the test function
#' \deqn{\mathbf{T} = (\overline{T}_1({\bf r}), \overline{T}_2({\bf r}), \dots , \overline{T}_J({\bf r}))}{T = (\bar{T}_1(r), \bar{T}_2(r), ..., \bar{T}_J(r))}
#' where \eqn{\overline{T}_i({\bf r})}{\bar{T}_i(r)} is a vector of mean values of functions in the group i.
#' This can be done by choosing the summaryfun \code{means}. This test assumes that the variances are
#' equal across the groups of functions. If the variances are not equal, one should use the summaryfun
#' \code{studmeans} instead. To account for unequal variances, a different test vector is considered.
#' From each of the \eqn{\overline{T}_i({\bf r})}{\bar{T}_i(r)} the average of all the test functions
#' without the test functions in the ith group is reduced and further this remainder is rescaled by
#' variances, see details in Mrkvička et al. (2016).
#'
#' An alternative is to test the equivalent hypothesis
#' \deqn{H_0 : \mu_i(r) - \mu_j(r) = 0, i=1,\dots,J-1, j=1,\dots,J.}{H0: \mu_i(r) - \mu_j(r) = 0, i=1,...,J-1, j=1,...,J.}
#' This test corresponds to the post-hoc test done usually after an ANOVA test is significant, but
#' it can be directed tested by mean of the combined rank test (Mrkvička et al., 2017), if the
#' test vector is taken to consist of the differences of the group averages of test functions, namely
#' \deqn{\mathbf{T'} = (\overline{T}_1({\bf r})-\overline{T}_2({\bf r}),
#' \overline{T}_1({\bf r})-\overline{T}_3({\bf r}), \dots , \overline{T}_{J-1}({\bf r})-\overline{T}_J({\bf r})).}{T' = (\bar{T}_1(r)-\bar{T}_2(r), \bar{T}_1(r)-\bar{T}_3(r), ..., \bar{T}_{J-1}(r)-\bar{T}_J(r)).}
#' The summaryfun option \code{contrasts} can be used to perform the test based on this test vector.
#' This again assumes that the variances are equal across the groups of functions. To deal with
#' unequal variances, the differences are rescaled, which corresponds to the summaryfun option
#' \code{studcontrasts}. See Mrkvička et al. (2016) for further details.
#'
#' @param nsim The number of random permutations.
#' @param x The original data (an array of functions). Typically a matrix or a data frame,
#' also \code{\link[fda.usc]{fdata}} objects allowed. Number of rows in x should correspond
#' to the number of groups, and each row should correspond to a function.
#' @param groups The original groups (a factor vector representing the assignment to groups).
#' @param summaryfun Possible values are "means", "studmeans", "contrasts", "studcontrasts".
#' See description for their meaning.
# Note: Possibly add a some arguments to specify which contrasts should be used.
# (Try to find our how this is usually done in R, in ordinary anova.)
#' @param alpha The significance level of the test.
#' @param n.aver If summaryfun is either "studmeans" or "studcontrasts", there is a possibility to
#' use variances smoothed by appying moving average to the estimated variance. n.iter determines
#' how many values on each side do contribute (incl. value itself).
#' @param mirror The complement of the argument circular of \code{\link[stats]{filter}}.
#' @param saveperm Logical. If TRUE, then the functions from permutations are saved to the attribute
#' simfuns.
#' @export
#' @references
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2016)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#'
#' Mrkvička, T., Myllymäki, M., and Hahn, U. (2017).
#' Multiple Monte Carlo testing, with applications in spatial point processes.
#' Statistics and Computing 27 (5): 1239-1255. doi:10.1007/s11222-016-9683-9
#' @examples
#' data(rimov)
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#'
#' nsim <- 2499 # for exploring
#'
#' res <- graph.fanova(nsim=nsim, x=rimov, groups=groups, summaryfun="means")
#' plot(res)
#' res2 <- graph.fanova(nsim=nsim, x=rimov, groups=groups, summaryfun="contrasts")
#' plot(res2)
graph.fanova <- function(nsim, x, groups, summaryfun, alpha=0.05, n.aver = 1L, mirror = FALSE, saveperm=FALSE) {
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(!(class(x) %in% c("matrix", "data.frame", "array", "fdata"))) stop("x is not a valid object.\n")
  if(nrow(x) != length(groups)) stop("The number of rows in x and the length of groups should be equal.\n")

  summaryfun <- spatstat::pickoption("sumf", summaryfun, c(means = "means",
                                                           mean = "means",
                                                           studmeans = "studmeans",
                                                           studmean = "studmeans",
                                                           contrasts = "contrasts",
                                                           cont = "contrasts",
                                                           studcontrasts = "studcontrasts",
                                                           studcont = "studcontrasts"))
  # setting that 'summaryfun' is a function
  switch(summaryfun, 
         means = {fun = means},
         studmeans = {fun = studmeans},
         contrasts = {fun = contrasts},
         studcontrasts = {fun = studcontrasts}
         )

  obs <- fun(x, groups, n.aver, mirror)
  # simulations by permuting to which groups the functions belong to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE), n.aver, mirror))
  # labels for comparisons (rownames(sim) is the same as names(obs))
  complabels <- unique(names(obs))

  cset <- create_curve_set(list(r = rep(1:ncol(x), times=length(complabels)),
                                obs = obs,
                                sim_m = sim))
  res_rank <- rank_envelope(cset, alpha=alpha, lexo=TRUE, alternative="two.sided")

  res <- structure(list(ranktest = res_rank,
                        summaryfun = summaryfun,
                        labels = complabels,
                        call = match.call()), class = "graph.fanova")
  if(saveperm) attr(res, "simfuns") <- sim
  res
}

#' Print method for the class 'graph.fanova'
#' @usage \method{print}{graph.fanova}(x, ...)
#'
#' @param x an 'graph.fanova' object
#' @param ... Ignored.
#'
#' @method print graph.fanova
#' @export
print.graph.fanova <- function(x, ...) {
  cat("Functional rank ANOVA with graphical interpretation (plot the object for it!).\n")
  print(x[['ranktest']])
}

#' Plot method for the class 'graph.fanova'
#' @usage \method{plot}{graph.fanova}(x, plot_style="ggplot2", separate_yaxes = TRUE, labels=x[['labels']], ...)
#'
#' @param x An 'graph.fanova' object
#' @param plot_style Either "basic" or "ggplot2".
#' @inheritParams env_ggplot
#' @param ... Additional parameters to be passed to \code{\link{plot.envelope_test}}.
#' @method plot graph.fanova
#' @export
plot.graph.fanova <- function(x, plot_style="ggplot2", separate_yaxes = TRUE, labels=x[['labels']], ...) {
  plot_style <- spatstat::pickoption("ptype", plot_style, c(basic = "basic",
                                                            b = "basic",
                                                            ggplot2 = "ggplot2",
                                                            ggplot = "ggplot2",
                                                            g = "ggplot2"))
  switch(plot_style,
         basic = {
           plot.envelope_test(x[['ranktest']], separate_yaxes=separate_yaxes, labels = NULL, ...) # ignoring labels
         },
         ggplot2 = {
           plot.envelope_test(x[['ranktest']], plot_style="ggplot2", separate_yaxes=separate_yaxes, labels = labels, ...)
         })
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
#' @param equalvar Logical. Default TRUE means that the traditional F-values are used.
#' In the case of FALSE, corrected F-values are used. The current implementation uses
#' \code{\link[stats]{lm}} to get the corrected F-values.
#' @export
#' @references
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2016)
#' A one-way ANOVA test for functional data with graphical interpretation.
#' arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)
#' @examples
#' \donttest{
#' data(rimov)
#' groups <- factor(c(rep(1, times=12), rep(2, times=12), rep(3, times=12)))
#' res <- frank.fanova(nsim=9999, x=rimov, groups=groups)
#' plot(res)
#' }
frank.fanova <- function(nsim, x, groups, alpha=0.05, equalvar=TRUE) {
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(!(class(x) %in% c("matrix", "data.frame", "array", "fdata"))) stop("x is not a valid object.\n")
  if(nrow(x) != length(groups)) stop("The number of rows in x and the length of groups should be equal.\n")

  if(equalvar) fun <- Fvalues
  else fun <- corrFvalues
  obs <- fun(x, groups)
  # simulations by permuting to which groups the functions belong to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE)))

  cset <- create_curve_set(list(r = 1:length(obs), obs = obs, sim_m = sim))
  # Perform the rank envelope test and return just the p-value
  res <- rank_envelope(cset, alpha=alpha, lexo=TRUE, alternative="greater")

  res
}
