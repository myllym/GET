context('p-values')



# Helpers.

case_no_ties <- function(position, n_sample) {
    list(argument=list(obs=position,
                       # Permuting the values for a harder test.
                       sim_vec=sample(seq_len(n_sample)[-position])),
         # Calculation modified from a separate source:
         # B. V. North, D. Curtis, and P. C. Sham, "A Note on the
         # Calculation of Empirical P Values from Monte Carlo Procedures,"
         # Am J Hum Genet, vol. 71, no. 2, pp. 439–441, 2002.
         p_reference=(n_sample - position + 1L) / n_sample)
}

do_no_ties_test <- function(position, n_sample) {
    func <- 'estimate_p_value.default'

    case_l <- case_no_ties(position, n_sample)
    argument_l <- case_l[['argument']]
    p_reference <- case_l[['p_reference']]

    p_midrank <- do.call(func, modifyList(argument_l,
                                          list(ties = 'midrank')))
    p_random <- do.call(func, modifyList(argument_l,
                                         list(ties = 'random')))
    p_conservative <- do.call(func, modifyList(argument_l,
                                               list(ties = 'conservative')))
    p_liberal <- do.call(func, modifyList(argument_l,
                                          list(ties = 'liberal')))

    expect_that(p_midrank, equals(p_reference))
    expect_that(p_random, equals(p_reference))
    expect_that(p_conservative, equals(p_reference))
    expect_that(p_liberal, equals(p_reference))
}



# Actual tests.

test_that('p-value is correct; no ties; one case', {
    expect_that(estimate_p_value(4, c(1, 2, 3)), equals(0.25))
})

test_that('p-value is correct; no ties; generated cases', {
    lapply(seq_len(10L), do_no_ties_test, 10L)
    NULL
})

test_that('p-value is correct; several ties; one case', {
    data_sample <- 2L
    mc_samples <- c(1L, 2L, 2L, 2L, 3L)

    p_reference_liberal <- 2 / 6
    p_reference_conservative <- 5 / 6
    p_reference_midrank <-
        1 - (rank(c(data_sample, mc_samples), ties.method='average')[1]
             - 1L) / 6
    p_reference_random_check <- function(x) {
        isTRUE(all.equal(x, 2 / 6)) ||
        isTRUE(all.equal(x, 3 / 6)) ||
        isTRUE(all.equal(x, 4 / 6)) ||
        isTRUE(all.equal(x, 5 / 6))
    }

    p_liberal <- estimate_p_value(data_sample, mc_samples, 'liberal')
    p_conservative <- estimate_p_value(data_sample, mc_samples,
                                       'conservative')
    p_midrank <- estimate_p_value(data_sample, mc_samples, 'midrank')
    p_random <- estimate_p_value(data_sample, mc_samples, 'random')

    expect_that(p_liberal, equals(p_reference_liberal))
    expect_that(p_conservative, equals(p_reference_conservative))
    expect_that(p_midrank, equals(p_reference_midrank))
    expect_that(p_reference_random_check(p_random), is_true())
})
