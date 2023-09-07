# Regress the given data (true or permuted) against the full model and
# get an effect of interest at all tau values in a matrix.
# @param Y True or permuted values of Y.
# @param df data
# nameinteresting = names of the interesting coefficients
#' @importFrom stats coef
genCoefmeans.rq <- function(Y, df, taus, formula.full, nameinteresting, ...) {
  anyNAY <- anyNA(Y)
  effect.a <- sapply(seq_along(taus), function(i) {
    if(is.vector(Y)) {
      df$response__ <- Y
    } else {
      df$response__ <- Y[,i]
    }
    if(anyNAY) {
      df <- df[!is.na(df$response__),]
    }
    M.full <- quantreg::rq(formula.full, data=df, tau=taus[i], model=FALSE, ...)
    allcoef <- coef(M.full)
    unlist(allcoef[nameinteresting])
  })
  if(length(nameinteresting)==1) {
    name <- names(effect.a)[1]
    dim(effect.a) <- c(length(effect.a), 1)
    dimnames(effect.a) = list(NULL, name)
    effect.a
  } else t(effect.a)
}


#' Graphical functional quantile regression
#'
#' Non-parametric graphical tests of significance in quantile regression
#'
#'
#' See the paper (Section 4) for details on the permutation methods.
#'
#'
#' @inheritParams graph.flm
#' @param data data.frame where to look for variables.
#' @param taus The quantiles to be used.
#' @param permutationstrategy The permutation strategy to be used. See details.
#' @param rq.args Additional arguments passed to \code{rq}.
#' @param contrasts Passed directly to \code{rq}.
#' @return A \code{global_envelope} or \code{combined_global_envelope} object,
#' which can be printed and plotted directly.
#' @export
#' @references
#' Mrkvička, T., Konstantinou, K., Kuronen, M. and Myllymäki, M. () Global quantile regression.
#'
#' Myllymäki, M and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#'
#' Freedman, D., & Lane, D. (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298. doi:10.2307/1391660
#' @importFrom parallel mclapply
#' @importFrom parallel parLapply
#' @importFrom parallel clusterEvalQ
#' @importFrom stats as.formula model.matrix update.formula
#' @examples
#' if(require("quantreg", quietly=TRUE)) {
#'   data("stackloss")
#'   res <- graph.rq(nsim=19,
#'     formula.full=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
#'     formula.reduced=stack.loss ~ Water.Temp,
#'     taus=seq(0.1, 0.9, length=10), permutationstrategy="remove quantile",
#'     data=stackloss, typeone="fwer")
#'   plot(res)
#' }
#'
graph.rq <- function(nsim, formula.full, formula.reduced, taus, typeone = c("fwer", "fdr"),
                     data = NULL, contrasts = NULL,
                     permutationstrategy = c("Freedman-Lane", "Freedman-Lane+remove zeros",
                                             "simple",
                                             "within nuisance",
                                             "remove location", "remove location scale",
                                             "remove quantile"),
                     savefuns = FALSE, rq.args = NULL, lm.args = NULL, GET.args = NULL,
                     mc.cores = 1L, mc.args = NULL, cl = NULL) {
  if(!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required, but not installed.")
  }
  typeone <- check_typeone(typeone, missing(typeone))

  if(missing(permutationstrategy)) stop("Permutation strategy has to be specified.")
  permutationstrategy = match.arg(permutationstrategy)

  check_isnested(formula.full, formula.reduced)
  if(nsim < 1) stop("Not a reasonable value of nsim.")

  genCoef <- genCoefmeans.rq
  ylab <- expression(italic(hat(beta)[i](tau)))

  nameinteresting <- setdiff(
    colnames(model.matrix(formula.full, data=data, contrasts.arg=contrasts)),
    colnames(model.matrix(formula.reduced, data=data, contrasts.arg=contrasts)))

  # Fit the full model to the data and obtain the coefficients
  # genCoef will use response__ as the response variable
  response_variable <- formula.full[[2]]
  formula.full <- update.formula(formula.full, response__ ~ .)
  Y <- eval(response_variable, data)

  # TODO: contrasts for variables not in reduced model will cause warning, but maybe are needed if there are interactions
  # Maybe make another contrasts where the extra variables are not included

  if(startsWith(permutationstrategy, "Freedman-Lane") || permutationstrategy=="remove quantile") {
    #-- Freedman-Lane procedure
    # Fit the reduced model at each argument value to get the fitted values and residuals
    fit <- do.call(quantreg::rq, c(list(formula.reduced, data=data, tau=taus, model=FALSE, contrasts=contrasts), rq.args))
    fitted.m <- fit$fitted.values
    res.m <- fit$residuals
    if(length(taus)==1) {
      fitted.m <- matrix(fitted.m, ncol=1)
      res.m <- matrix(res.m, ncol=1)
    }

    if(permutationstrategy=="remove quantile") {
      # Note that formula.full is actually formula.interesting
      Y <- res.m
      fitted.m <- 0
      formula.full <- update.formula(formula.full, as.formula(paste("~.-(", as.character(formula.reduced)[3], ")")))
    } else if(permutationstrategy=="Freedman-Lane+remove zeros") {
      # Find the indices of the zero residuals that will be removed
      # Since there can be arbitrarily many zeros, take a random sample
      names.reduced <- if(length(taus)==1) {names(coef(fit))} else {rownames(coef(fit))}
      zerostoremove <- length(names.reduced) - 1
      if(zerostoremove >= 1) {
        zerores <- abs(res.m) < 1e-10
        stopifnot(all(colSums(zerores) >= zerostoremove))
        # Set residuals to be removed to NA, any rows with NA are removed in gencoef
        for(i in seq_along(taus)) {
          ind <- sample(which(zerores[,i]), zerostoremove, replace=FALSE)
          res.m[ind, i] <- NA
        }
      }
    }
    fit <- NULL
  } else if(permutationstrategy=="simple") {
    fitted.m <- 0
    res.m <- matrix(Y, ncol=1)
  } else if(permutationstrategy=="within nuisance") {
    Nuisance_name = as.character(formula.reduced[[3]])
    if(length(Nuisance_name) != 1) stop("formula_reduced may only contain 1 variable if 'within nuisance' permutation strategy is used.")
    Nuisance.levels <- unique(data[[Nuisance_name]])
    ## TODO: Maybe check that the nuisance has a reasonable amount of data on each level?

    # permutation strategy ( Works only for one categorical covariate )
    Yperm <- function() {
      result <- Y_temp <- Y
      for( i in Nuisance.levels){
        Index_i = which(data[[Nuisance_name]]==i)
        permutation <- Index_i[sample.int(length(Index_i))]
        result[Index_i] = Y_temp[permutation]
      }
      result
    }
  } else if(startsWith(permutationstrategy, "remove location")) {
    # These are implemented by computing the corresponding residuals and
    # then performing the permutation as simple permutation for the residuals and
    # the formula where nuisance are removed.

    #simulated data
    fit <- do.call(lm, c(list(formula.reduced, data = data, contrasts=contrasts),lm.args)) # Fit the reduced model

    res.m <- fit$residuals #residuals

    #NP + permutation
    if(permutationstrategy == "remove location scale") {
      data.temp = data
      data.temp$abs_residuals = abs(res.m)
      formula.res <- update.formula(formula.reduced, abs_residuals ~ .)
      fit.abs<- do.call(lm, c(list( formula.res,data = data.temp, contrasts=contrasts),lm.args))
      fitted.m = fit.abs$fitted.values
      res.m <- res.m / fitted.m
    }

    formula.full <- update.formula(formula.full, as.formula(paste("~.-(", as.character(formula.reduced)[3], ")")))
    Y <- res.m
    res.m <- matrix(res.m, ncol=1)
    fitted.m <- 0
  }

  obs <- do.call(genCoef, c(list(Y, data, taus, formula.full, nameinteresting, contrasts=contrasts), rq.args))

  if(permutationstrategy != "within nuisance") {
    # Simulations by permuting the residuals + calculate the coefficients
    Yperm <- function() { # Permutation
      permutation <- sample(1:nrow(res.m), size=nrow(res.m), replace=FALSE)
      # Permute the residuals (rows in res.m) and create new 'y'
      fitted.m + res.m[permutation, ]
    }
  }
  loopfun2 <- function(i) {
    # Regress the permuted data against the full model and get a new effect of interest
    do.call(genCoef, c(list(Yperm(), data, taus, formula.full, nameinteresting, contrasts=contrasts), rq.args))
  }
  if(is.null(cl)) sim <- do.call(parallel::mclapply, c(list(X=1:nsim, FUN=loopfun2, mc.cores=mc.cores), mc.args))
  else sim <- parLapply(cl, 1:nsim, loopfun2)
  dn <- dimnames(sim[[1]])
  # sim <- simplify2array(sim, except=NULL) # except is only available starting from 4.2.0
  sim <- simplify2array(sim)
  if(is.null(dimnames(sim))) {
    dim(sim) <- c(1,1,length(sim))
    dimnames(sim) <- c(dn, list(NULL))
  }
  complabels <- colnames(obs)

  csets <- NULL
  for(i in 1:ncol(obs)) {
    simi <- array(sim[,i,], dim=dim(sim)[c(1,3)]) # Extract the slice even if some dimensions are 1
    csets[[complabels[i]]] <- create_curve_set(list(r = taus,
                                                    obs = obs[,i],
                                                    sim_m = simi))
  }
  switch(typeone,
         fwer = {
           res <- do.call(global_envelope_test,
                          c(list(curve_sets=csets, alternative="two.sided", nstep=1), GET.args))
         },
         fdr = {
           res <- do.call(fdr_envelope,
                          c(list(curve_sets=csets, alternative="two.sided"), GET.args))
         })
  attr(res, "method") <- "Graphical Quantile Regression" # Change method name
  attr(res, "permutationstrategy") <- permutationstrategy
  attr(res, "labels") <- complabels
  # Re-define the default ylab
  res <- envelope_set_labs(res, ylab=ylab)
  attr(res, "call") <- match.call()
  if(savefuns) attr(res, "simfuns") <- csets
  res
}
