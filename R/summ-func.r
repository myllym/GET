#' Helper function for both \code{\link{estimate_K_f}} and
#' \code{\link{estimate_L_f}}.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_unmarked Whether to include the unmarked K or L
#'   in the result.
#' @param do_besags_L Whether to calculate the K_f or the L_f function.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'unmarked'
#'   contains the unmarked K or L function, if asked for.
call_summ_func <- function(pattern, do_besags_L,
                           edge_correction = 'translate', mtf_name = 'm',
                           r_max = NULL, r_vec = NULL, calc_unmarked = TRUE,
                           use_biased_lambda2 = FALSE) {
    got_req <- require(marksummary)
    if (!got_req) {
        stop('marksummary must be installed for call_summ_func.')
    }

    if (length(mtf_name) != 1L || mtf_name %in% '1') {
        stop('mtf_name must be a single string and not "1".')
    }
    if (length(calc_unmarked) != 1L || !is.logical(calc_unmarked) ||
        !is.finite(calc_unmarked)) {
        stop('calc_unmarked must be TRUE or FALSE.')
    }

    if (calc_unmarked) {
        call_mtf_name <- c('1', mtf_name)
    } else {
        call_mtf_name <- mtf_name
    }
    summ <- summ_func(pattern = pattern, edge_correction = edge_correction,
                     mtf_name = call_mtf_name, r_max = r_max, r_vec = r_vec,
                     do_besags_L = do_besags_L,
                     use_biased_lambda2 = use_biased_lambda2)

    if (do_besags_L) {
        name_prefix <- 'L'
    } else {
        name_prefix <- 'K'
    }

    a <- summ[['a']]
    obs_name <- paste(name_prefix, mtf_name, sep = '_')
    res <- list(r = summ[['r']], obs = a[obs_name, , drop = TRUE])

    if (calc_unmarked) {
        unmarked_name <- paste(name_prefix, '1', sep = '_')
        res[['unmarked']] <- a[unmarked_name, , drop = TRUE]
    }

    res
}

#' Estimate K_f.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_unmarked Whether to include the unmarked K in the result.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'unmarked'
#'   contains the unmarked K function, if asked for.
#' @export
estimate_K_f <- function(pattern, mtf_name = 'm', r_max = NULL,
                         r_vec = NULL, calc_unmarked = TRUE,
                         edge_correction = 'translate',
                         use_biased_lambda2 = FALSE) {
    do_besags_L <- FALSE
    res <- call_summ_func(pattern = pattern, do_besags_L = do_besags_L,
                          edge_correction = edge_correction,
                          mtf_name = mtf_name, r_max = r_max,
                          r_vec = r_vec, calc_unmarked = calc_unmarked,
                          use_biased_lambda2 = use_biased_lambda2)
}

#' Estimate L_f.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_unmarked Whether to include the unmarked L in the result.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'unmarked'
#'   contains the unmarked L function, if asked for.
#' @export
estimate_L_f <- function(pattern, mtf_name = 'm', r_max = NULL,
                         r_vec = NULL, calc_unmarked = TRUE,
                         edge_correction = 'translate',
                         use_biased_lambda2 = FALSE) {
    do_besags_L <- TRUE
    res <- call_summ_func(pattern = pattern, do_besags_L = do_besags_L,
                          edge_correction = edge_correction,
                          mtf_name = mtf_name, r_max = r_max,
                          r_vec = r_vec, calc_unmarked = calc_unmarked,
                          use_biased_lambda2 = use_biased_lambda2)
}


#' Helper function for both \code{\link{random_labelling_K_f}} and
#' \code{\link{random_labelling_L_f}}.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_theo Whether to calculate the theoretical summary function
#'   as well.
#' @param do_besags_L Whether to calculate the K_f or the L_f function.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'unmarked'
#'   contains the unmarked K or L function, if asked for.
call_random_labelling <- function(pattern, do_besags_L, mtf_name = 'm',
                                  r_max = NULL, r_vec = NULL, n_sim = 999L,
                                  calc_theo = TRUE,
                                  edge_correction = 'translate',
                                  method = 'permute',
                                  use_biased_lambda2 = FALSE) {
    got_req <- require(marksummary)
    if (!got_req) {
        stop('marksummary must be installed for call_random_labelling.')
    }

    if (length(mtf_name) != 1L || mtf_name %in% '1') {
        stop('mtf_name must be a single string and not "1".')
    }
    if (length(calc_theo) != 1L || !is.logical(calc_theo) ||
        !is.finite(calc_theo)) {
        stop('calc_theo must be TRUE or FALSE.')
    }

    if (calc_theo) {
        call_mtf_name <- c('1', mtf_name)
    } else {
        call_mtf_name <- mtf_name
    }
    summ <- summ_func_random_labelling(pattern,
                                       edge_correction = edge_correction,
                                       mtf_name = call_mtf_name,
                                       n_perm = n_sim, r_max = r_max,
                                       r_vec = r_vec,
                                       do_besags_L = do_besags_L,
                                       method = method,
                                       use_biased_lambda2 =
                                           use_biased_lambda2)

    if (do_besags_L) {
        name_prefix <- 'L'
    } else {
        name_prefix <- 'K'
    }

    a <- summ[['a']]
    obs_name <- paste(name_prefix, mtf_name, sep = '_')
    sim_idx <- !(dimnames(a)[['orig_and_perm']] %in% 'original')
    res <- list(r = summ[['r']],
                obs = a['original', obs_name, , drop = TRUE],
                sim_m = t(a[sim_idx, obs_name, , drop = TRUE]))

    if (calc_theo) {
        unmarked_name <- paste(name_prefix, '1', sep = '_')
        res[['theo']] <- a['original', unmarked_name, , drop = TRUE]
    }

    res[['is_residual']] <- FALSE

    res <- create_curve_set(res)
    res
}

#' Generate K_f curves for a random labelling test.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_theo Whether to calculate the theoretical summary function
#'   as well.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'theo'
#'   contains the unmarked K function, if asked for.
#' @export
random_labelling_K_f <- function(pattern, mtf_name = 'm', r_max = NULL,
                                 r_vec = NULL, n_sim = 999L,
                                 calc_theo = TRUE,
                                 edge_correction = 'translate',
                                 method = 'permute',
                                 use_biased_lambda2 = FALSE) {
    do_besags_L <- FALSE
    res <- call_random_labelling(pattern, do_besags_L = do_besags_L,
                                 mtf_name = mtf_name, r_max = r_max,
                                 r_vec = r_vec, n_sim = n_sim,
                                 calc_theo = calc_theo,
                                 edge_correction = edge_correction,
                                 method = method,
                                 use_biased_lambda2 = use_biased_lambda2)
}

#' Generate L_f curves for a random labelling test.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @param calc_theo Whether to calculate the theoretical summary function
#'   as well.
#' @return A list with either two or three components. 'obs' has the summary
#'   function that was asked for. 'r' contains the radius values. 'theo'
#'   contains the unmarked L function, if asked for.
#' @export
random_labelling_L_f <- function(pattern, mtf_name = 'm', r_max = NULL,
                                 r_vec = NULL, n_sim = 999L,
                                 calc_theo = TRUE,
                                 edge_correction = 'translate',
                                 method = 'permute',
                                 use_biased_lambda2 = FALSE) {
    do_besags_L <- TRUE
    res <- call_random_labelling(pattern, do_besags_L = do_besags_L,
                                 mtf_name = mtf_name, r_max = r_max,
                                 r_vec = r_vec, n_sim = n_sim,
                                 calc_theo = calc_theo,
                                 edge_correction = edge_correction,
                                 method = method,
                                 use_biased_lambda2 = use_biased_lambda2)
}
