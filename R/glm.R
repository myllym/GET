# Preliminary checks for the graph.fglm and frank.fglm
fglm.checks <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL) {
  # Preliminary checks
  vars <- all.vars(formula.full)
  vars.reduced <- all.vars(formula.reduced)
  if(!all(vars.reduced %in% vars)) stop("The reduced model includes something extra, not in the full model.\n")
  if(nsim < 1) stop("Not a reasonable value of nsim.\n")
  if(class(curve_sets)[1] != "list") {
    curve_sets <- list(Y=curve_sets)
    if(vars[1] != "Y") stop("The formula should be off the form Y ~ .... where Y is the response.\n")
  }
  if(!( all(vars %in% names(curve_sets)) | all(vars[-1] %in% names(factors)) )) stop("The variables in the formula not found in the given data (curve_sets and factors).\n")
  if(!all(lapply(curve_sets, class) %in% c("curve_set", "fdata"))) stop("The components of curve_sets do not have a valid class.\n")
  curve_sets <- lapply(curve_sets, convert_fdata)
  if(!all(lapply(curve_sets[['obs']], is.matrix))) stop("The curve_set must include data functions (sim_m ignored).\n")
  curve_sets <- check_curve_set_dimensions(curve_sets)
  # Put the functions and factors into data.l
  data.l <- list()
  data.l[[1]] <- t(curve_sets[[vars[1]]][['obs']]) # -> each row corresponds to a data function
  # The argument values
  r <- curve_sets[[vars[1]]][['r']]
  Nfunc <- nrow(data.l[[1]]) # Number of functions
  nr <- ncol(data.l[[1]]) # Number of argument values
  if(length(curve_sets) > 1) { # Factors provided in the curve_sets
    for(i in 2:length(vars)) data.l[[i]] <- t(curve_sets[[vars[i]]][['obs']])
  }
  else {
    if(is.null(factors)) stop("No factors provided.\n")
    else {
      if(class(factors)[1] != "data.frame") stop("Invalid factors argument.\n")
      if(nrow(factors) != nrow(data.l[[1]])) stop("The dimensions of Y and factors do not match.\n")
      # Expand the factors to each argument value
      for(i in 2:length(vars)) data.l[[i]] <- matrix(factors[,vars[i]], nrow=nrow(factors), ncol=nr)
    }
  }
  names(data.l) <- vars
  list(data.l=data.l, r=r, Nfunc=Nfunc, nr=nr)
}

# Regress the given data (true or permuted) against the full model and
# get an effect of interest at all r values as a long vector.
# data.l = a list containing data (Y and factors), all variables in formula.full
# nameinteresting = names of the interesting coefficients
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefmeans <- function(data.l, formula.full, nameinteresting, ...) {
  nr <- ncol(data.l[[1]])
  effect.a <- matrix(0, nrow=nr, ncol=length(nameinteresting))
  dimnames(effect.a) <- list(NULL, nameinteresting)
  for(i in 1:nr) {
    df <- as.data.frame(lapply(data.l, FUN = function(x) x[,i]))
    M.full <- stats::lm(formula.full, data=df, ...)
    allcoef <- unlist(stats::dummy.coef(M.full))
    effect.a[i,] <- allcoef[nameinteresting]
  }
  res <- as.vector(effect.a) # to a vector
  names(res) <- rep(nameinteresting, each=nr)
  res
}

# The constrasts of the interesting effects
#' @importFrom stats lm
#' @importFrom stats dummy.coef
genCoefcontrasts <- function(data.l, formula.full, nameinteresting, ...) {
  nr <- ncol(data.l[[1]])
  k <- length(nameinteresting)
  effect.a <- matrix(0, nrow=nr, ncol=k)
  dimnames(effect.a) <- list(NULL, nameinteresting)
  for(i in 1:nr) {
    df <- as.data.frame(lapply(data.l, FUN = function(x) x[,i]))
    M.full <- stats::lm(formula.full, data=df, ...)
    allcoef <- unlist(stats::dummy.coef(M.full))
    effect.a[i,] <- allcoef[nameinteresting]
  }
  # contrasts
  cont <- NULL
  for(i in 1:(k-1)) for(j in (i+1):k)  {
    ct <- effect.a[,i] - effect.a[,j] # coef_i - coef_j
    names(ct) <- rep(paste(nameinteresting[i], nameinteresting[j], sep="-"), nr)
    cont <- c(cont, ct)
  }
  cont
}

# General F-values from lm-model
#' @importFrom stats lm
#' @importFrom stats anova
genFvalues <- function(data.l, formula.full, formula.reduced, ...) {
  nr <- ncol(data.l[[1]])
  Fvalues <- vector(length=nr)
  for(i in 1:nr) {
    df <- as.data.frame(lapply(data.l, FUN = function(x) x[,i]))
    M.full <- stats::lm(formula = formula.full, data = df, ...)
    M.reduced <- stats::lm(formula = formula.reduced, data = df, ...)
    Anova.res <- stats::anova(M.reduced, M.full)
    Fvalues[i] <- Anova.res$F[2]
  }
  Fvalues
}


#' Graphical functional GLM
#'
#' Graphical functional general linear model (GLM)
#'
#'
#' The function \code{graph.fglm} performs the graphical functional GLM of Mrkvička et al. (2019).
#' This is a nonparametric graphical test of significance of a covariate in functional GLM.
#' The test is able to find not only if the factor of interest is significant, but also which
#' functional domain is responsible for the potential rejection.
#' In the case of functional multi-way main effect ANOVA or functional main effect ANCOVA models,
#' the test is able to find which groups differ (and where they differ).
#' In the case of functional factorial ANOVA or functional factorial ANCOVA models,
#' the test is able to find which combination of levels (which interactions) differ (and where they differ).
#' The described tests are global envelope tests applied in the context of GLMs.
#' The Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
#' (to obtain the simulations under the null hypothesis of "no effects");
#' consequently, the test approximately achieves the desired significance level.
#'
#' @inheritParams graph.fanova
#' @inheritParams dg.global_envelope_test
#' @param formula.full The formula specifying the general linear model,
#' see \code{formula} in \code{\link[stats]{lm}}.
#' @param formula.reduced The formula of the reduced model with nuisance factors only.
#' @param curve_sets A named list of sets of curves giving the dependent variable (Y), and
#' possibly additionally all the factors. The dimensions of the elements should
#' match with each other, i.e. the factor values should be given for each argument value
#' and each function.
#' @param factors A data frame of factors. An alternative way to specify factors when they
#' are constant for all argument values. The number of rows of the data frame should be equal
#' to the number of curves. Each column should specify the values of a factor.
#' @param ... Additional arguments to be passed to \code{\link[stats]{lm}}.
#' @param GET.args A named list of additional arguments to be passed to \code{\link{global_envelope_test}}.
#' @param mc.args A named list of additional arguments to be passed to \code{\link[parallel]{mclapply}}.
#' Only relevant if \code{mc.cores} is more than 1.
#' @return A \code{global_envelope} object.
#' @export
#' @references
#' Mrkvička, T., Roskovec, T. and Rost, M. (2019). A nonparametric graphical tests of significance in functional GLM. arXiv:1902.04926 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @importFrom stats lm
#' @importFrom stats predict.lm
#' @importFrom stats dummy.coef
#' @importFrom parallel mclapply
#' @examples
#' data(GDPtax)
#' factors.df <- data.frame(Group = GDPtax$Group, Tax = GDPtax$Profittax)
#' res.tax_within_group <- graph.fglm(nsim = 999,
#'                                   formula.full = Y~Group+Tax+Group:Tax,
#'                                   formula.reduced = Y~Group+Tax,
#'                                   curve_sets = GDPtax$curve_set,
#'                                   factors = factors.df)
#' plot(res.tax_within_group, plot_style="ggplot2")
#'
graph.fglm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                       summaryfun = c("means", "contrasts"),
                       savefuns = FALSE, ..., GET.args = NULL, mc.cores = 1, mc.args = NULL) {
  # Set up the contrasts
  options(contrasts = c("contr.sum", "contr.poly"))
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- fglm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors)
  summaryfun <- match.arg(summaryfun)
  # setting that 'fun' is a function
  switch(summaryfun, 
         means = {genCoef = genCoefmeans},
         contrasts = {genCoef = genCoefcontrasts}
  )

  # Fit the full model at the first argument value to get the variables
  df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,1]))
  mod.full <- lm(formula.full, data=df, ...)
  namecoef.full <- names(unlist(stats::dummy.coef(mod.full)))
  # Freedman-Lane procedure
  # Fit the reduced model at each argument value
  fitted.m <- res.m <- matrix(0, X$Nfunc, X$nr)
  for(i in 1:X$nr) {
    df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,i])) # create the data.frame at the ith argument value
    mod.red <- stats::lm(formula.reduced, data=df, ...)
    # Save predictions and residuals
    fitted.m[,i] <- mod.red$fitted.values
    res.m[,i] <- mod.red$residuals
  }
  # Names of the coefficients in the reduced model
  namecoef.red <- names(unlist(stats::dummy.coef(mod.red)))
  # Change the names of the coefficients in the reduced model,
  # if the full model includes discrete factors, but the reduced model not
  if(length(mod.red$xlevels) == 0 & length(mod.full$xlevels) > 0) {
    namecoef.red <- paste(namecoef.red, ".", namecoef.red, sep="")
  }
  # The interesting coefficients
  nameinteresting <- setdiff(namecoef.full, namecoef.red)
  cat("The inspected coefficients are: \n")
  cat(nameinteresting)
  # Fit the full model to the data and obtain the effects
  obs <- genCoef(X$data.l, formula.full, nameinteresting, ...)
  # Simulations by permuting the residuals
  loopfun <- function(i, ...) {
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    X$data.l[[1]] <- fitted.m + res.m[permutation, ]
    # Regress the permuted data against the full model and get a new effect of interest
    genCoef(X$data.l, formula.full, nameinteresting, ...)
  }
  sim <- do.call(parallel::mclapply, c(list(X=1:nsim, FUN=loopfun, mc.cores=mc.cores), mc.args))
  sim <- sapply(sim, function(x) x, simplify="array")
  complabels <- unique(names(obs))

  cset <- create_curve_set(list(r = rep(X$r, times=length(complabels)),
                                obs = obs,
                                sim_m = sim))
  res <- do.call(global_envelope_test, c(list(curve_sets=cset, alternative="two.sided"), GET.args))
  attr(res, "method") <- "Graphical functional GLM" # Change method name
  attr(res, "labels") <- complabels
  attr(res, "call") <- match.call()
  if(savefuns) attr(res, "simfuns") <- sim
  res
}

#' F rank functional GLM
#'
#' Multiple testing in permutation inference for the general linear model (GLM)
#'
#'
#' The function \code{frank.fglm} performs
#' a nonparametric test of significance of a covariate in the functional GLM.
#' Similarly as in the graphical functional GLM (\code{\link{graph.fglm}}),
#' the Freedman-Lane algorithm (Freedman and Lane, 1983) is applied to permute the functions
#' (to obtain the simulations under the null hypothesis of "no effects");
#' consequently, the test approximately achieves the desired significance level.
#' In contrast to the graphical functional GLM, the F rank functional GLM is based on the F-statistics,
#' that are calculated at each argument value of the functions.
#' The global envelope test is applied to the observed and simulated F-statistics.
#' The test is able to find if the factor of interest is significant and also which
#' argument values of the functional domain are responsible for the potential rejection.
#'
#' @inheritParams graph.fglm
# Freedman-Lane procedure (Freedman and Lane, 1983, p. 385)
frank.fglm <- function(nsim, formula.full, formula.reduced, curve_sets, factors = NULL,
                       ..., GET.args = NULL, mc.cores = 1, mc.args = NULL) {
  # Preliminary checks and formulation of the data to suitable form for further processing
  X <- fglm.checks(nsim=nsim, formula.full=formula.full, formula.reduced=formula.reduced,
                   curve_sets=curve_sets, factors=factors)

  # Freedman-Lane procedure
  # Fit the reduced model at each argument value
  fitted.m <- res.m <- matrix(0, X$Nfunc, X$nr)
  for(i in 1:X$nr) {
    df <- as.data.frame(lapply(X$data.l, FUN = function(x) x[,i])) # create the data.frame at the ith argument value
    mod.red <- stats::lm(formula.reduced, data=df, ...)
    # Save predictions and residuals
    fitted.m[,i] <- mod.red$fitted.values
    res.m[,i] <- mod.red$residuals
  }
  # Calculate the F-statistic for the data
  obs <- genFvalues(X$data.l, formula.full, formula.reduced, ...)
  # Simulations by permuting the residuals + calculate F-values
  loopfun <- function(i, ...) {
    permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
    # Permute the residuals (rows in res.m) and create new 'y'
    X$data.l[[1]] <- fitted.m + res.m[permutation, ]
    # Regress the permuted data against the full model and get a new effect of interest
    genFvalues(X$data.l, formula.full, formula.reduced, ...)
  }
  sim <- do.call(parallel::mclapply, c(list(X=1:nsim, FUN=loopfun, mc.cores=mc.cores), mc.args))
  sim <- sapply(sim, function(x) x, simplify="array")

  cset <- create_curve_set(list(r = X$r, obs = obs, sim_m = sim))
  do.call(global_envelope_test, c(list(curve_sets=cset, alternative="greater"), GET.args))
}
