context('p-value correctness compared to devtests')

test_that('spruces gives the same result', {
    pattern <- spatstat::spruces

    # Random seed.
    seed <- 1234L

    mtf_name <- 'm'
    # Number of permutations.
    n_sim <- 999L
    r_min <- NULL
    r_max <- NULL
    r_vec <- NULL
    use_L <- TRUE
    edge_correction <- 'translate'

    # Change both at the same time.
    measure <- 'max'
    devtests_measure <- 'max'

    # Change both at the same time.
    scaling <- 'env'
    devtests_scaling <- 'Enve'

    # Do not change these.
    use_theo <- TRUE
    method <- 'permute'
    use_biased_lambda2 <- FALSE
    ties <- 'random'



    if (use_L) {
        base_letter <- 'L'
    } else {
        base_letter <- 'K'
    }
    transf <- paste(base_letter, devtests_scaling, sep = '')
    devtests_row <- paste(transf, '_', mtf_name, sep = '')
    transform_l <- devtests::PrepareTransformations(transf)

    set.seed(seed)
    p_all_devtests <- devtests:::RunPermutationTest(
                          pattern, n.perm = n_sim,
                          transform.l = transform_l, r.min = r_min,
                          r.max = r_max, use.theoretical = use_theo,
                          edge.correction = edge_correction)
    p_devtests <- p_all_devtests[devtests_row, devtests_measure]



    set.seed(seed)
    p_spptest <- random_labelling_test(pattern, mtf_name = mtf_name,
                                       n_sim = n_sim, r_min = r_min,
                                       r_max = r_max, r_vec = r_vec,
                                       measure = measure, scaling = scaling,
                                       use_L = use_L,
                                       edge_correction = edge_correction,
                                       use_theo = use_theo, method = method,
                                       use_biased_lambda2 =
                                           use_biased_lambda2,
                                       ties = ties)



    expect_that(p_spptest, equals(p_devtests))
    expect_that(p_spptest, is_identical_to(p_devtests))
})
