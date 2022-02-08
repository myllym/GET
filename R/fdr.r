# Returns the last component of a vector
vlast <- function(x) {
  x[length(x)]
}

# Helper functions for FDR envelopes

# Calculate Qi's for each resample for the Yekutieli & Benjamini (1999) FDR estimation:
# If u gives \rho_u^0(p) and pi0=1, then the Q values correspond to Yekutieli & Benjamini (1999) 'u'-version.
# If u gives m*p, then it is our modification - which we finally did not study.
#
# The expectation EQ is the mean of the given Qs.
# @param u Either the level of the envelope, nr*gamma, or the beta-quantile
# @param R0s Rejections for the functions from the null for the gamma envelope
# @param Robs The number of rejections for the observed function for the gamma envelope
# @param nr Number of hypotheses
Qcalc <- function(u, R0s, Robs, pi0) {
  num <- pi0 * R0s
  dem <- pi0 * (R0s - u) + Robs
  Q <- num/dem
  Q[R0s == 0] <- 0
  Q[R0s > 0 & Robs < pi0*u] <- 1 #Q[dem <= 0 & R0s != 0] <- 1
  Q
}
Qcalc_m <- function(u_vect, rejs_sim, rejs_obs, pi0) {
  nsim2 <- nrow(rejs_sim)
  u_mat <- matrix(rep(u_vect, times=nsim2), byrow=TRUE, nrow=nsim2)
  rejs_obs_mat <- matrix(rep(rejs_obs, times=nsim2), byrow=TRUE, nrow=nsim2)
  num <- pi0 * rejs_sim
  dem <- pi0 * (rejs_sim - u_mat) + rejs_obs_mat
  Q <- num/dem
  Q[rejs_sim == 0] <- 0
  Q[rejs_sim > 0 & rejs_obs_mat < pi0*u_mat] <- 1
  Q
}

# Original Yekutieli & Benjamini (1999) Q's
# Eq. (6) of Mrkvicka and Myllymäki (2022)
# @param u1, u2 The levels of the envelope, u1 = nr*gamma, u2 = beta-quantile
Qcalc_YB <- function(u1, u2, R0s, Robs, pi0) {
  num <- pi0 * R0s
  dem <- pi0 * (R0s - u1) + Robs
  Q <- num/dem
  Q[R0s == 0] <- 0
  Q[R0s > 0 & Robs-u2 < u1] <- 1
  Q
}
Qcalc_m_YB <- function(u1_vect, u2_vect, rejs_sim, rejs_obs, pi0) {
  nsim2 <- nrow(rejs_sim)
  u1_mat <- matrix(rep(u1_vect, times=nsim2), byrow=TRUE, nrow=nsim2)
  u2_mat <- matrix(rep(u2_vect, times=nsim2), byrow=TRUE, nrow=nsim2)
  rejs_obs_mat <- matrix(rep(rejs_obs, times=nsim2), byrow=TRUE, nrow=nsim2)
  num <- pi0 * rejs_sim
  dem <- pi0 * (rejs_sim - u1_mat) + rejs_obs_mat
  Q <- num/dem
  Q[rejs_sim == 0] <- 0
  Q[rejs_sim > 0 & rejs_obs_mat-u2_mat < pi0*u1_mat] <- 1
  Q
}


# Calculate E(R0) and R.obs, and E(Q) if curve_set2 given, for the rank envelopes (gamma = 1,2,3,4,...)
# @param curve_set2 If curve_set2 given, then also E(Q) and E(R0) are calculated from this set. If not,
# then E(R0) is calculated using the theoretical expectation.
# @param fast TRUE means that different ranking of curves are calculated through other rankings.
#   FALSE calculates thoroughly, and is a version to compare to.
fdr_rejections_rank <- function(curve_set, alternative, curve_set2=NULL, pi0=1, beta=0.05, fast=TRUE) {
  curve_set <- convert_envelope(curve_set)

  obs_curve <- curve_set_1obs(curve_set)
  sim_curves <- data_or_sim_curves(curve_set) # Exclude the data function
  #all_curves <- data_and_sim_curves(curve_set) # Include the data function
  nsim <- curve_set_nfunc(curve_set) - 1
  nr <- curve_set_narg(curve_set) # m in the fdr paper

  # Count the numbers of rejections R0_j for data and for all simulations j,
  # j=1,...,nsim2 and for all envelopes l=1,2,3,...,floor(nsim/2)
  if(alternative=="two.sided") {
    max_l <- floor(nsim/2)
    mult <- 2
  }
  else {
    max_l <- floor(3/4*nsim) # In the one-sided case, go up to 25% envelope).
    mult <- 1
  }

  if(!is.null(curve_set2)) {
    curve_set2 <- convert_envelope(curve_set2)

    # From the second set of simulations, calculate R0 (the numbers of rejections under the null)
    sim_curves2 <- data_or_sim_curves(curve_set2) # Exclude the data function (if it exists)
    # Note: It has been checked that the curves curve_set and curve_set2 have corresponding dimension, and obs_curve is the same
    all_curves2 <- data_and_sim_curves(curve_set2) # The data function and the second set of simulated functions
    nsim2 <- curve_set_nfunc(curve_set2) - 1

    # Q_u calculation
    # rejs_sim = matrix of rejections (rows~resamples, columns~rank envelopes up to max_l)
    # rejs_obs = vector of rejections (~rank envelopes)
    Q0_all_calc <- function(rejs_sim, rejs_obs, beta) {
      u1 <- nr*mult*(1:max_l)/nsim # m * tilde(gamma)
      u2 <- apply(rejs_sim, MARGIN = 2, FUN=quantile, probs=1-beta) # Parametric version: R_sd <- apply(rejs_sim, MARGIN = 2, FUN=sd); u2 <- qnorm(beta, mean=u1, sd=R_sd)
      Qm_YB <- Qcalc_m_YB(u1, u2, rejs_sim, rejs_obs, pi0)
      Qm_YBu <- Qcalc_m(u2, rejs_sim, rejs_obs, pi0)
      list(Q_YB=apply(Qm_YB, MARGIN = 2, FUN=mean),
           Q_YBu=apply(Qm_YBu, MARGIN = 2, FUN=mean))
    }

    if(fast) {
      all_curves <- rbind(obs_curve, sim_curves2, sim_curves)

      # Calculate pointwise ranks for each argument value (r)
      calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
      for(i in 1:nr) {
        sim_curves[,i] <- calc_pointwiserank(sim_curves[,i]) # overwriting curves by their ranks
        all_curves[,i] <- calc_pointwiserank(all_curves[,i])
        all_curves2[,i] <- calc_pointwiserank(all_curves2[,i])
      }
      # Ranks of the simulated curves among themselves from which the envelopes are constructed
      simranks <- sim_curves
      # Ranks of the data curve and all simulated curves from curve_set2 among all the curves
      datasim2ranks <- all_curves[1:(1+nsim2),]
      # Ranks of the data curve and all simulated curves from curve_set2 among curve_set curves
      datasim2ranks <- datasim2ranks - (all_curves2 - 1)

      rejs <- matrix(nrow=1+nsim2, ncol=max_l)
      #Q <- matrix(nrow=nsim2, ncol=max_l)
      for(i in 1:max_l) {
        rejected <- datasim2ranks <= i
        rejs[,i] <- rowSums(rejected)
        #Q[,i] <- Qcalc(nr*mult*i/nsim, rejs[-1,i], rejs[1,i], pi0)
      }
      R0 <- apply(rejs[-1,], MARGIN = 2, FUN=mean)
      R.obs <- rejs[1,]
      #Q0 <- apply(Q, MARGIN = 2, FUN=mean)
      #Q0_u <- Q0_u_calc(rejs[-1,], rejs[1,], beta)
      Q0all <- Q0_all_calc(rejs[-1,], rejs[1,], beta)
      dataranks <- datasim2ranks[1,] # Save the dataranks for qvalue calculations
    }
    else {
      # Calculate pointwise ranks for each argument value (r) for the functions in the first set of curves
      calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
      simranks <- sim_curves
      for(i in 1:nr) {
        simranks[,i] <- calc_pointwiserank(sim_curves[,i]) # NOT overwriting curves by their ranks
      }

      # Calculate ranks for data curve
      all_curves <- rbind(obs_curve, sim_curves)
      for(i in 1:nr) {
        all_curves[,i] <- calc_pointwiserank(all_curves[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
      }
      dataranks <- all_curves[1,]

      # and all curves from the second set of simulations
      sim2_ranks <- sim_curves2
      for(j in 1:nsim2) {
        all_curves2 <- rbind(sim_curves2[j,], sim_curves)
        for(i in 1:nr) {
          all_curves2[,i] <- calc_pointwiserank(all_curves2[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
        }
        sim2_ranks[j,] <- all_curves2[1,]
      }

      # Count the numbers of rejections
      # for data
      R.obs <- sapply(1:max_l, FUN = function(i) sum(dataranks <= i))
      # for the second set of simulations
      rejs <- matrix(nrow=nsim2, ncol=max_l)
      for(i in 1:max_l) {
        rejected <- sim2_ranks <= i
        rejs[,i] <- rowSums(rejected)
      }
      R0 <- apply(rejs, MARGIN = 2, FUN=mean)
      Q0all <- Q0_all_calc(rejs, R.obs, beta)
    }
  }
  else { # The case of just one set of curves (curve_set2 = NULL)
    # Use the theoretical expectation for ER0:
    all_curves <- data_and_sim_curves(curve_set) # Include the data function

    # Calculate pointwise ranks for each argument value (r)
    calc_pointwiserank <- find_calc_pointwiserank("rank", alternative)
    for(i in 1:nr) {
      all_curves[,i] <- calc_pointwiserank(all_curves[,i]) # overwriting curves by their ranks. FIXME(for computational reasons)? Only the rank of data needed here.
    }
    dataranks <- all_curves[1,]
    R.obs <- sapply(1:max_l, FUN = function(i) sum(dataranks <= i))

    # The expected numbers of rejections R0 for all envelopes l=1,2,3...,max_l
    R0 <- nr * mult * (1:max_l)/nsim
    Q0all <- NULL
  }

  list(R0=R0, R.obs=R.obs, Q0=Q0all, dataranks=dataranks)
}

# Storey & Tibshirani (2001)
fdr_step13 <- function(R.obs, R0) {
  R0 / pmax(R.obs, 1)
}

# Find the index from rejs (R.obs, R0) for which the fdr is below the given limit
find_ind <- function(R.obs, R0, limit, leftright="right") {
  switch(leftright,
    right = { # The last one that is still within the limit
      max(0, vlast(which(fdr_step13(R.obs, R0) <= limit)))
    },
    left = { # The last one before fdr first goes above the limit
      max(0, (which(fdr_step13(R.obs, R0) > limit))[1]) - 1
    }
  )
}
find_ind_YB <- function(Q, limit, leftright="right") {
  switch(leftright,
         right = { # The last one that is still within the limit
           max(0, vlast(which(Q <= limit)))
         },
         left = { # The last one before fdr first goes above the limit
           max(0, (which(Q > limit))[1]) - 1
         }
  )
}

# Calculate the lower and upper bound of the envelope
fdr_envelope_rank <- function(ind, curve_set) {
  # ind corresponds to the l of the lth minimal and maximal pointwise values
  sim_curves <- data_or_sim_curves(curve_set)
  nr <- curve_set_narg(curve_set)
  nsim <- dim(sim_curves)[1]
  LB <- UB <- array(0, nr)
  for(i in 1:nr){
    Hod <- sort(sim_curves[,i])
    LB[i] <- Hod[ind]
    UB[i] <- Hod[nsim-ind+1]
  }
  list(LB=LB, UB=UB)
}

# Given the index of the 'l'ind'th rank envelope, and the number of data rejections for that ind_Robs (to handle a special case),
# find the gamma between suitable the indices [ind, ind+1].
# If curve_set2 == NULL, calculates the rejections for the FDRest choice "ST".
# If curve_set2 is provided, then also for "YB" and "YBu" (which require arguments pi0 and beta).
# Note that anyway only the one for which 'ind' is given is useful and the function using this function
# should take care of that.
fdr_rejections_between_two_rank <- function(ind, ind_Robs, limit, curve_set, curve_set2=NULL,
                                            alternative="two.sided", pi0=NULL, beta=0.05) {
  failed <- FALSE
  obs_curve <- curve_set_1obs(curve_set)
  nr <- curve_set_narg(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1
  if(!is.null(curve_set2)) {
    sim_curves2 <- data_or_sim_curves(curve_set2) # The second set of simulations; calculate the rejections for these
    nsim2 <- dim(sim_curves2)[1]
  }
  noutside <- function(curve, lower_bound, upper_bound, alternative) {
    switch(alternative,
      two.sided = { sum(curve < lower_bound | curve > upper_bound) },
      less = { sum(curve < lower_bound) },
      greater = { sum(curve > upper_bound) }
    )
  }
  if(alternative == "two.sided") mult <- 2 else mult <- 1
  if(ind == 0) { # The value of gamma is between 0 and 1
    # Find the envelope of rank 1; THE envelope is wider than this
    e2 <- fdr_envelope_rank(1, curve_set)
    # e2l <- e2u <- array(0, nr)
    # for(i in 1:nr) {
    #   Hod <- sort(sim_curves[,i])
    #   e2l[i] <- Hod[ind+1]
    #   e2u[i] <- Hod[nsim-ind]
    # }
    # Define the gamma values
    # Smallest gamma needs to be smaller than nsim * limit / ( nr * mult )
    d <- max(0.0001, min(0.0075, nsim * limit / ( nr * mult )))
    gammas <- seq(0+d/2, 1-d/2, by=d)
    #gammas <- seq(0, 1, length=130); gammas <- gammas[-c(1, length(gammas))]
    # Calculate rejections for each gamma
    e <- list()
    R.obs_new <- R0_new <- R0_new_sd <- Q_new <- Qu_new <- vector(length=length(gammas))
    for(i in 1:length(gammas)) {
      if(alternative != "greater") LB <- e2$LB + log(gammas[i]) * (e2$UB - e2$LB) else LB <- -Inf
      if(alternative != "less") UB <- e2$UB - log(gammas[i]) * (e2$UB - e2$LB) else UB <- Inf
      R.obs_new[i] <- noutside(obs_curve, LB, UB, alternative)
      if(!is.null(curve_set2)) {
        R0s <- sapply(1:nsim2, FUN = function(j) { noutside(sim_curves2[j,], LB, UB, alternative) })
        R0_new[i] <- mean(R0s)
        #R0_new_sd[i] <- sd(R0s)
        # Q-calculation for the Yekutieli & Benjamini (1999) FDR estimation
        u1 <- nr*mult*gammas[i]/nsim
        u2 <- quantile(R0s, probs=1-beta) # qnorm(1-beta, mean=u1, sd=R0_new_sd[i])
        Q_new[i] <- mean(Qcalc_YB(u1, u2, R0s, R.obs_new[i], pi0))
        Qu_new[i] <- mean(Qcalc(u2, R0s, R.obs_new[i], pi0))
      }
      else
        R0_new[i] <- nr*mult*gammas[i]/nsim
      e[[i]] <- list(LB=LB, UB=UB)
    }
  }
  else { # The value of gamma lies between ind and ind + 1
    # Calculate the envelopes between which THE envelope lies
    sim_curves <- data_or_sim_curves(curve_set)
    e1l <- e1u <- e2l <- e2u <- array(0, nr)
    for(i in 1:nr) {
      Hod <- sort(sim_curves[,i])
      e1l[i] <- Hod[ind]
      e1u[i] <- Hod[nsim-ind+1]
      e2l[i] <- Hod[ind+1]
      e2u[i] <- Hod[nsim-ind]
    }
    switch(alternative,
           "two.sided" = {},
           "less" = { e1u <- e2u <- Inf },
           "greater" = { e1l <- e2l <- -Inf })
    # Define gamma values
    gammas <- seq(ind, ind+1, length=18)
    gammas <- gammas[-c(1, length(gammas))]
    # Calculate rejections for each gamma
    e <- list()
    #rejs_new <- matrix(nrow=nsim2, ncol=length(gammas))
    R.obs_new <- R0_new <- R0_new_sd <- Q_new <- Qu_new <- vector(length=length(gammas))
    for(i in 1:length(gammas)) {
      if(alternative != "greater") LB <- e1l + (gammas[i] - ind)*(e2l-e1l) else LB <- -Inf
      if(alternative != "less") UB <- e1u - (gammas[i] - ind)*(e1u-e2u) else UB <- Inf
      R.obs_new[i] <- noutside(obs_curve, LB, UB, alternative)
      if(!is.null(curve_set2)) {
        R0s <- sapply(1:nsim2, FUN = function(j) { noutside(sim_curves2[j,], LB, UB, alternative) })
        #rejs_new[,i] <- R0s
        R0_new[i] <- mean(R0s)
        R0_new_sd[i] <- sd(R0s)
        # Q-calculation for the Yekutieli & Benjamini (1999) FDR estimation
        u1 <- nr*mult*gammas[i]/nsim
        u2 <- quantile(R0s, probs=1-beta) # qnorm(beta, mean=u1, sd=R0_new_sd[i])
        Q_new[i] <- mean(Qcalc_YB(u1, u2, R0s, R.obs_new[i], pi0))
        Qu_new[i] <- mean(Qcalc(u2, R0s, R.obs_new[i], pi0))
      }
      else
        R0_new[i] <- nr*mult*gammas[i]/nsim
      e[[i]] <- list(LB=LB, UB=UB)
    }
  }
  #list(gammas=gammas, R0=R0_new, R.obs=R.obs_new, Q0=Q_new, Q0_u=Qu_new)
  #-- Find the gamma
  find_gamma <- function(ind_new) {
    if(ind_new == 0) {
      if(ind != 0) {
        gamma <- ind
        e <- list(LB=e1l, UB=e1u)
        r1 <- ind_Robs # Robs[ind]
      }
      else {
        warning("Envelope was not found. Returning the largest envelope.")
        gamma <- gammas[1]
        e <- e[[1]]
        r1 <- R.obs_new[1]
        failed <- TRUE
      }
    }
    else {
      gamma <- gammas[ind_new]
      e <- e[[ind_new]]
      r1 <- R.obs_new[ind_new]
    }
    list(gamma=gamma, e=e, r1=r1, failed=failed)
  }
  # Indices for the three methods
  ind_new_ST <- find_ind(R.obs_new, R0_new, limit)
  if(!is.null(curve_set2)) {
    ind_new_YB <- find_ind_YB(Q_new, limit)
    ind_new_YBu <- find_ind_YB(Qu_new, limit)
    resYB <- find_gamma(ind_new_YB)
    resYBu <- find_gamma(ind_new_YBu)
  }
  else resYB <- resYBu <- NULL

  list(ST=find_gamma(ind_new_ST), YB=resYB, YBu=resYBu)
}


general_FDRenvelope_algorithm <- function(algorithm = c("ATSE", "IATSE"),
                                          FDRest= c("ST", "YB", "YBu"),
                                          alpha, alternative, curve_set, curve_set2 = NULL, beta = 0.05,
                                          rejs1=NULL) {
  pi0 <- 1
  gammastar <- step2 <- NULL
  nr <- curve_set_narg(curve_set)
  nsim <- curve_set_nfunc(curve_set) - 1

  switch(algorithm,
         ATSE = {
           limit1 <- function(alpha) alpha/(1+alpha)
           limit2 <- function(alpha, pi0) alpha/pi0
           get_pi0est <- function(r1, gamma=NULL) { # r1 (gamma ignored) from step1
             (nr-r1)/nr
           }
         },
         IATSE = {
           limit1 <- function(alpha) alpha
           limit2 <- function(alpha, pi0) alpha/vlast(pi0)
           get_pi0est <- function(r1, gamma) { # r1 & gamma from step1
             if(alternative == "two.sided") gammastar <- 2*gamma/nsim
             else gammastar <- gamma/nsim
             pi0[1] <- min(1, (nr-r1+nr*gammastar)/nr)
             pi0[2] <- min(1, (nr-r1+pi0[1]*nr*gammastar)/nr)
             j <- 2
             while(pi0[j-1]-pi0[j] > 10^-10) {
               j <- j+1
               pi0[j] <- min(1, (nr-r1+pi0[j-1]*nr*gammastar)/nr)
             }
             pi0
           }
         })
  switch(FDRest,
         ST = {
           find_ind_local <- function(rejs, limit) find_ind(rejs$R.obs, rejs$R0, limit)
           if(!is.null(curve_set2)) warning("Using the second set of simulations for the E(R0) estimation instead of the theoretical formula.")
         },
         YB = {
           find_ind_local <- function(rejs, limit) find_ind_YB(rejs$Q0$Q_YB, limit)
           if(is.null(curve_set2)) warning("The second set of simulations not provided for the E(Q) estimation, using an approximation/ratio of expectations.")
         },
         YBu = {
           find_ind_local <- function(rejs, limit) find_ind_YB(rejs$Q0$Q_YBu, limit)
           if(is.null(curve_set2)) warning("The second set of simulations not provided for the E(Q) estimation, using an approximation/ratio of expectations.")
         }
  )

  # Calculate E(R0), R.obs, E(Q)
  if(is.null(rejs1)) rejs1 <- fdr_rejections_rank(curve_set, alternative, curve_set2, pi0=1, beta=beta) # curve_set NULL (ER0) or given (also EQ)
  # Find the largest rejection region G* for which E(Q) <= alpha/(1+alpha)
  ind1 <- find_ind_local(rejs1, limit1(alpha))
  step1 <- fdr_rejections_between_two_rank(ind1, rejs1$R.obs[ind1], limit1(alpha), curve_set, curve_set2, alternative,
                                           pi0=pi0, beta=beta) # pi0=1
  # Number of rejected hypotheses
  r1 <- step1[[FDRest]]$r1
  if(!(r1==0 | r1==nr) | (algorithm == "IATSE" & r1==nr)) {
    # Estimate \pi_0
    pi0 <- get_pi0est(r1, step1[[FDRest]]$gamma)
    # Find the largest rejection region G for which E(R0) / max(R.obs, 1) <= alpha/pi0
    if(FDRest != "ST") { # Need to recalculate the Q's again, because they depend on pi0
      rejs1 <- fdr_rejections_rank(curve_set, alternative, curve_set2, pi0=vlast(pi0), beta=beta) # rejs2, but using the same object name
    }
    ind2 <- find_ind_local(rejs1, limit2(alpha, pi0))
    step2 <- fdr_rejections_between_two_rank(ind2, rejs1$R.obs[ind2], limit2(alpha, pi0), curve_set, curve_set2, alternative,
                                             pi0=vlast(pi0), beta=beta)
  }
  else # Take the first step envelope as the envelope
    step2 <- step1

  list(step2=step2[[FDRest]], pi0=pi0, rejs=rejs1, step1=step1[[FDRest]]) # what else? Note: rejs includes dataranks
}


# Calculates the fdr envelope for one curve_set
individual_fdr_envelope <- function(curve_set, curve_set2 = NULL, alpha = 0.05,
                                    alternative = c("two.sided", "less", "greater"),
                                    algorithm=c("IATSE", "ATSE"), FDRest = c("ST", "YB", "YBu"),
                                    beta = 0.05, savefdr=FALSE) {
  isenvelope <- inherits(curve_set, "envelope")
  alternative <- match.arg(alternative)
  algorithm <- match.arg(algorithm)
  FDRest <- match.arg(FDRest)
  if(!is.numeric(alpha) || (alpha < 0 | alpha >= 1)) stop("Unreasonable value of alpha.")
  if(!is.numeric(beta) || (beta <= 0 | beta >= 1)) stop("Unreasonable value of beta.")
  picked_attr <- pick_attributes(curve_set, alternative=alternative) # Saving for attributes / plotting purposes

  curve_set <- convert_envelope(curve_set)
  if(!curve_set_is1obs(curve_set)) stop("The curve_set must contain an observed function.")
  obs_curve <- curve_set_1obs(curve_set)

  if(!is.null(curve_set2)) {
    curve_set2 <- convert_envelope(curve_set2)
    if(curve_set_is1obs(curve_set2))
      if(!identical(obs_curve, curve_set_1obs(curve_set2)))
        stop("The observed functions in the two sets of curves differ.")
    if(!identical(curve_set_narg(curve_set), curve_set_narg(curve_set2)))
      stop("The number of hypotheses in the two sets of curves differ.")
  }

  # Go through the algorithm to find the 'indices' for constructing envelopes
  inds <- general_FDRenvelope_algorithm(algorithm=algorithm, FDRest=FDRest,
                                alpha, alternative, curve_set, curve_set2=NULL, beta=beta)

  res <- make_envelope_object(type=algorithm, curve_set, LB=inds$step2$e$LB, UB=inds$step2$e$UB,
                              T_0=rowMeans(curve_set[['funcs']][,-1]),
                              picked_attr=picked_attr, isenvelope=isenvelope,
                              Malpha=NULL, alpha=alpha, distance=NULL)
  # Change the "method" attribute
  attr(res, "method") <- "FDR envelope"

  class(res) <- c("fdr_envelope", class(res))
  attr(res, "rej") <- obs_curve < inds$step2$e$LB | obs_curve > inds$step2$e$UB
  attr(res, "gamma") <- inds$step2$gamma
  attr(res, "pi0") <- inds$pi0
  if(savefdr) {
    attr(res, "R0") <- inds$rejs$R0
    attr(res, "R") <- inds$rejs$R.obs
    attr(res, "Q0") <- inds$rejs$Q0
    attr(res, "Q0u") <- inds$rejs$Q0u
  }
  res
}

#' The FDR envelope
#'
#' Calculate the FDR envelope based on ATSE or IATSE algorithm.
#'
#' @param curve_sets A \code{curve_set} (see \code{\link{create_curve_set}}) or an
#'   \code{envelope} object of \pkg{spatstat} containing the observed function
#'   and the functions from which the envelope is to be constructed.
#'   Alternatively, a list of appropriate objects can be given.
#' @param curve_sets2 A \code{curve_set} giving another set of simulations from which
#' the numbers of rejections under the null are calculated. (Alternatively, a list similarly
#' as \code{curve_sets}.) If not provided, then expectation under H0 is used instead following
#' Storey (2002) and Storey and Tibshirani (2003).
#' @param type The type of the global envelope with current options for 'rank', 'erl', 'cont', 'area'.
#' See \code{\link{forder}} for description of the measures.
#' @inheritParams global_envelope_test
#' @param method Either "ATSE" or "IATSE" standing for adaptive two stage
#'   global envelope, and iteratively adaptive two stage global envelope methods,
#'   respectively.
#' @export
#' @examples
#' # Just for fun
#' if(require("spatstat", quietly=TRUE)) {
#'   # CSR test
#'   #=========
#'   X <- unmark(spruces)
#'   \donttest{nsim <- 2000 # Number of simulations}
#'   \dontshow{nsim <- 20}
#'   env1 <- envelope(X, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   transform = expression(.-r), # centering
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   res <- fdr_envelope(env1)
#'   plot(res)
#' }
#'
#' # A GLM example
#' data(rimov)
#' \donttest{nsim <- 1000 # Number of simulations}
#' \dontshow{nsim <- 20}
#' sims1 <- graph.flm(nsim=nsim,
#'                  formula.full = Y~Year,
#'                  formula.reduced = Y~1,
#'                  curve_sets = list(Y=rimov), factors = data.frame(Year = 1979:2014),
#'                  savefuns = TRUE)
#' res <- fdr_envelope(attr(sims1, "simfuns")[[1]])
#' plot(res)
#'
fdr_envelope <- function(curve_sets, alpha = 0.05, curve_sets2=NULL,
                         alternative = c("two.sided", "less", "greater"),
                         algorithm = c("IATSE", "ATSE"), FDRest = c("ST", "YB", "YBu"),
                         beta = 0.05) {
  # Consider the case of several curves
  if(length(class(curve_sets)) == 1 && class(curve_sets) == "list") {
    if(length(curve_sets) > 1) { # Combine
      nr <- lapply(curve_sets, curve_set_narg)
      nfuns <- length(nr)
      is2d <- curve_set_is2d(curve_sets[[1]])
      csetnames <- names(curve_sets)
      curve_sets <- combine_curve_sets(curve_sets, equalr=FALSE) # Allow also different numbers of hypotheses
      if(!is.null(curve_sets2)) {
        if(!identical(nr, lapply(curve_sets2, curve_set_narg)))
          stop("The numbers of hypotheses in the two sets of curves differ.")
        curve_sets2 <- combine_curve_sets(curve_sets2, equalr=FALSE)
      }
      # Test based on the combined sets
      res <- individual_fdr_envelope(curve_sets, curve_set2=curve_sets2, alpha=alpha,
                                     alternative=alternative,
                                     algorithm=algorithm, FDRest=FDRest,
                                     beta=beta, savefdr=FALSE)
      # Transform the envelope to a combined envelope
      idx <- unlist(lapply(1:nfuns, FUN = function(i) { rep(i, times=nr[i]) } ))
      # Split the envelopes to the original groups
      res_ls <- split(res, f=idx)
      # Redefine attributes method and rej
      rejs <- split(attr(res, "rej"), f=idx)
      for(i in 1:nfuns) {
        attr(res_ls[[i]], "method") <- paste0("1/", nfuns, "th of a combined FDR envelope")
        attr(res_ls[[i]], "rej") <- rejs[[i]]
      }
      # Set unreasonable attributes of individuals sets of curves to NULL
      anames <- c("alpha", "R0", "R", "pi0")
      anames <- anames[anames %in% names(attributes(res_ls[[1]]))]
      for(name in anames) {
        for(i in 1:nfuns) attr(res_ls[[i]], name) <- NULL
      }
      mostattributes(res_ls) <- attributes(res)
      attr(res_ls, "row.names") <- NULL
      if(!is.null(names(curve_sets))) names(res_ls) <- csetnames
      attr(res_ls, "method") <- "Combined FDR envelope"
      attr(res_ls, "nstep") <- 1
      class(res_ls) <- c("combined_fdr_envelope", "combined_global_envelope", "list")
      if(is2d)
        class(res_ls) <- c("combined_global_envelope2d", class(res_ls))
      return(res_ls)
    }
    else if(length(curve_sets) == 1) {
      curve_sets <- curve_sets[[1]] # unlist, and call individual_...
      curve_sets2 <- curve_sets[[1]]
    }
    else
      stop("The given list of curve_sets is empty.")
  }
  else {
    return(individual_fdr_envelope(curve_sets, curve_set2=curve_sets2, alpha=alpha,
                                   alternative=alternative,
                                   algorithm=algorithm, FDRest=FDRest,
                                   beta=beta, savefdr=FALSE))
  }
}



#' Print method for the class 'fdr_envelope'
#'
#' @param x An 'fdr_envelope' object
#' @param ... Ignored.
#' @export
print.fdr_envelope <- function(x, ...) {
  cat(attr(x, "method"), "based on", attr(x, "einfo")$nsim, "simulations:\n",
      "Number of rejected hypotheses:", sum(attr(x, "rej")), "\n",
      "Number of accepted hypotheses:", length(x[['r']]) - sum(attr(x, "rej")), "\n",
      "Total number of hypotheses:", length(x[['r']]), "\n")
}
