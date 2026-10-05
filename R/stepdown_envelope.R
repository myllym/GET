# Determine ordering to be used for stepdown envelope
find_most_extreme_erl <- function(ranks) {
  if(length(ranks) == 1) return(1)
  inds <- which(ranks == min(ranks))
  if(length(inds) == 1) return(inds)
  bestlen <- 1
  first <- 1
  len <- 1
  for(i in 2:length(inds)) {
    if(inds[i] == inds[i-1] + 1) {
      len <- len + 1
    } else {
      if(len > bestlen) {
        bestlen <- len
        first <- i - len
      }
      len <- 1
    }
  }
  if(len > bestlen) {
    bestlen <- len
    first <- i - len + 1
  }
  inds[first] + (bestlen - 1) %/% 2
}

find_most_extreme_area <- function(ranks) {
  R <- min(ceiling(ranks))
  if(length(ranks) == 1) return(1)
  inds <- which(ranks <= R)
  if(length(inds) == 1) return(inds)
  bestlen <- -1
  first <- 1
  len <- R - ranks[inds[1]]
  for(i in 2:length(inds)) {
    if(inds[i] == inds[i-1] + 1) {
      len <- len + (R - ranks[inds[i]])
    } else {
      if(len > bestlen) {
        bestlen <- len
        mid <- inds[first] - 1 + which.min(ranks[inds[first]:inds[i-1]])
      }
      first <- i
      len <- (R - ranks[inds[i]])
    }
  }
  if(len > bestlen) {
    bestlen <- len
    mid <- inds[first] - 1 + which.min(ranks[inds[first]:inds[i]])
  }
  mid
}

find_most_extreme <- function(ranks, measure) {
  switch (measure,
          erl = find_most_extreme_erl(ranks),
          area = find_most_extreme_area(ranks),
          stop("Unkown measure")
  )
}

determine_ordering <- function(ranks, measure) {
  if(measure %in% c("rank", "cont")) order(ranks)
  else {
    n <- length(ranks)
    ii <- 1:n
    o <- integer(n)
    for(i in 1:n) {
      oi <- find_most_extreme(ranks, measure)
      ranks <- ranks[-oi]
      o[i] <- ii[oi]
      ii <- ii[-oi]
    }
    o
  }
}

# Calculate pointwise ranks for the observed function
calc_ranks1 <- function(curve_set, measure = "rank", alternative) {
  all_curves <- data_and_sim_curves(curve_set)
  nr <- curve_set_narg(curve_set)
  # Calculate pointwise ranks for each argument value (r)
  calc_pointwiserank <- find_calc_pointwiserank(measure, alternative)
  ranks <- numeric(nr)
  for(i in 1:nr) {
    ranks[i] <- calc_pointwiserank(all_curves[,i])[1]
  }
  ranks
}

# Simple reference implementation that could be used in tests
internal_stepdown_envelope_slow <- function(curve_set, measure, alpha, alternative) {
  if(measure=="qdir") {
    probs <- c(alpha/2, 1-alpha/2)
    curve_set_o <- curve_set
    if(alternative %in% c("greater", "less")) {
      curve_set_o[['theo']] <- 0
    }
    curve_set_o <- residual(curve_set_o, use_theo = TRUE)
    curve_set_o <- scale_curves(curve_set_o, scaling = "qdir", probs = probs)
    ordering <- order(abs(curve_set_1obs(curve_set_o)), decreasing = TRUE)
  } else {
    ranks <- calc_ranks1(curve_set, measure, alternative)[1,]
    ordering <- determine_ordering(ranks, measure)
  }
  env <- NULL
  narg <- curve_set_narg(curve_set)
  inds <- 1:narg
  stage <- 1
  rejected <- logical(narg)
  while(stage <= narg){
    #cat(stage, " ")
    i <- ordering[stage]
    ii <- which(inds==i)
    inds <- inds[-ii]
    get1 <- global_envelope_test(curve_set, type=measure, alpha = alpha, alternative=alternative, ties="conservative")
    pvalue <- attr(get1, "p")
    if(pvalue > alpha + 1e-10) {
      stage <- narg
      env1 <- as.data.frame(get1)
      env1$rejected <- FALSE
    } else {
      env1 <- as.data.frame(get1)[ii,]
      env1$rejected <- TRUE
    }
    env1$pvalue <- pvalue
    env1$order <- stage
    env1$M_alpha <- attr(get1, "M_alpha")

    env <- rbind(env, env1)
    stage <- stage + 1
    curve_set <- curve_set[-ii,]
  }
  if(measure != "qdir") stopifnot(all.equal((env$obs > env$hi | env$obs < env$lo), env$rejected))
  env[order(env$r),]
}


# allpadj==TRUE computes all adjusted p-values
# allpadj==FALSE computes only adjusted p-values <= alpha and the correct envelope.
internal_stepdown_envelope <- function(curve_set, measure, envtype, alpha, alternative, erlhistn=6, debug=FALSE, allpadj=TRUE) {
  nr <- curve_set_narg(curve_set)
  Nfunc <- curve_set_nfunc(curve_set)

  if(startsWith(measure, "qdir") || measure=="unscaled") {
    probs <- c(alpha/2, 1-alpha/2)
    curve_set_res <- curve_set
    if(alternative %in% c("greater", "less")) {
      curve_set_res[['theo']] <- 0
    }
    T_0 <- get_T_0(curve_set, central="mean")

    if(measure == "unscaled" && alternative %in% c("greater", "less")) {
      T_0 <- rep(0, times=length(T_0))
    }

    curve_set_res <- residual(curve_set_res, use_theo = TRUE)
    curve_set_qdir <- scale_curves(curve_set_res, scaling = if(measure=="unscaled") "none" else "qdir", probs = probs)
    ordering <- order(abs(curve_set_1obs(curve_set_qdir)), decreasing = TRUE)

    all_curves <- unname(data_and_sim_curves(curve_set))
    all_curves_res <- unname(data_and_sim_curves(curve_set_res))
    allranks <- -abs(unname(data_and_sim_curves(curve_set_qdir))) # min to max
    measure <- "rank"
  } else {
    all_curves <- unname(data_and_sim_curves(curve_set))
    # Calculate pointwise ranks for each argument value (r)
    allranks <- all_curves
    calc_pointwiserank <- find_calc_pointwiserank(measure, alternative)
    for(i in 1:nr) {
      allranks[,i] <- calc_pointwiserank(all_curves[,i])
    }
    ordering <- determine_ordering(allranks[1,], measure)
  }

  orig_ranks <- allranks[1,]

  # Calculate measures from the pointwise ranks
  applyfuncs <- function(f) {
    sapply(1:nrow(allranks), function(i) { f(allranks[i,]) })
  }
  # Compute highest ranks (r) and their counts (c) for the n highest ranks
  # return r_1,-c_1,r_2,-c_2,r_3,-c_3,...,r_n,-c_n
  # ordering by the return value gives the same order as sorting by the ranks themselves
  # when there are no ties in the return values
  erlhist <- function(x) {
    y <- rle(sort(x))
    c(matrix(c(y$values[1:erlhistn], -y$lengths[1:erlhistn]), byrow=TRUE, nrow=2))
  }
  partialarea <- function(cont) {
    rank <- ceiling(min(cont))
    area <- sum(rank - cont[cont <= rank])
    c(rank, area)
  }
  index <- 1:nr
  padj <- numeric(nr)
  padj0 <- numeric(nr)
  ord0 <- integer(nr)
  ord_i <- 1
  LBounds <- UBounds <- numeric(nr)
  stage <- 1

  if(measure == "erl") {
    partial_forder <- if(erlhistn > 0) applyfuncs(erlhist)
    else applyfuncs(sort)
    droppedranks <- NULL
    partial_forder0 <- partial_forder
  }

  updatingenvelope <- TRUE
  while(ncol(allranks) > 0) { # nr
    if(measure %in% c("rank", "cont", "qdir", "unscaled")) {
      partial_forder <- applyfuncs(min)
    } else if(measure == "area") {
      partial_forder <- applyfuncs(partialarea)
    } else if(measure == "erl" && debug) {
      if(erlhistn > 0) partial_forder <- applyfuncs(erlhist)
      else partial_forder <- applyfuncs(sort)
    } else if(measure == "erl" && length(droppedranks)) {
      if(erlhistn > 0) {
        for(i in 1:nrow(allranks)) { # Nfunc
          j <- match(droppedranks[i,], partial_forder[c(TRUE, FALSE), i], nomatch=0) # 6
          if(any(j)) {
            if(partial_forder[2*j, i] < -1) {
              partial_forder[2*j, i] <- partial_forder[2*j, i] + 1
            } else {
              partial_forder[,i] <- erlhist(allranks[i,]) # nr log(nr)
            }
          }
          #stopifnot(all.equal(partial_forder[,i], erlhist(allranks[i,])))
        }
      } else {
        for(i in 1:nrow(allranks)) { # Nfunc
          dr <- droppedranks[i,]
          j <- match(dr, partial_forder0[, i], nomatch=0) # nr
          #stopifnot(!anyDuplicated(j))
          partial_forder[,i] <- partial_forder0[-j,i]
          #stopifnot(all.equal(partial_forder[,i], sort(as.vector(allranks[i,]))))
        }
      }
    }
    distance <- distance_from_partial_forder(measure, partial_forder, nr)
    i <- ordering[stage]
    ii <- which(index==i)

    droppedranks <- allranks[,ii, drop=FALSE]
    allranks <- allranks[,-ii, drop=FALSE]
    if(measure == "erl" && erlhistn == 0) {
      partial_forder0 <- partial_forder
      partial_forder <- partial_forder[-ii,, drop=FALSE]
    }

    Malpha <- critical(distance, alpha, Nfunc, TRUE)

    u <- -distance
    #pi <- estimate_p_value(u[1], u[-1], ties="conservative")
    pi <- (sum(u[-1] >= u[1]) + 1) / Nfunc
    padj[index[ii]] <- max(padj, pi)
    padj0[index[ii]] <- pi
    ord0[index[ii]] <- ord_i
    ord_i <- ord_i + 1
    stage <- stage + 1
    indexenv <- index[ii]
    index <- index[-ii]
    if(pi > alpha + 1e-15) {
      if(!allpadj)
        padj[index] <- 2
      indexenv <- c(index, indexenv)
    }
    if(updatingenvelope) {
      if(envtype=="functional") {
        j <- distance >= Malpha
        for(i in indexenv){
          lu <- range(all_curves[j,i])
          LBounds[i] <- lu[1]
          UBounds[i] <- lu[2]
        }
      } else if(envtype=="qdir") {
        for(i in indexenv) {
          quant_m <- quantile(all_curves_res[,i], probs)
          LBounds[i] <- T_0[i] - (-Malpha)*abs(quant_m[1])
          UBounds[i] <- T_0[i] + (-Malpha)*abs(quant_m[2])
        }
      } else if(envtype=="unscaled") {
        for(i in indexenv) {
          LBounds[i] <- T_0[i] - (-Malpha)
          UBounds[i] <- T_0[i] + (-Malpha)
        }
      } else if(envtype=="rank") {
        for(i in indexenv) {
          Hod <- sort(all_curves[,i])
          LBounds[i]<- Hod[Malpha]
          UBounds[i]<- Hod[Nfunc-Malpha+1]
        }
      } else {stop("envtype")}
    }
    if(pi > alpha + 1e-15) {
      updatingenvelope <- FALSE
      if(!allpadj)
        break
    }
  }
  rej <- padj <= alpha + 1e-15
  hi <- UBounds
  lo <- LBounds
  obs <- all_curves[1,]
  if(debug && envtype=="functional" && !isTRUE(all.equal((obs > hi | obs < lo), rej)))
    warning("IGI failed")
  # data.frame(lo, obs, hi, rej, padj, orig_ranks, ord0) |> mutate(ok=(obs > hi | obs < lo) == rej)
  # IGI does not hold if the number simulations is not "exact"
  list(padj=padj, lo=lo, hi=hi, obs=all_curves[1,], padj=padj, padj_not_correct=padj0, ord=ord0, rejected=rej)
}
