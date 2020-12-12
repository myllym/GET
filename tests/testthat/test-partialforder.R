test_that("correspondence_partial_forder_vs_forder", {
  data("abide_9002_23")
  cset <- frank.flm(nsim=19, formula.full = Y ~ Group + Sex + Age,
                    formula.reduced = Y ~ Group + Sex,
                    curve_sets = list(Y = abide_9002_23$curve_set),
                    factors = abide_9002_23$factors, savefuns = "return")
  p1 <- partial_forder(cset[1:100,], measure="area")
  p2 <- partial_forder(cset[-(1:100),], measure="area")
  expect_equal(combine_forder(list(p1, p2)), forder(cset, measure="area"))
  p1 <- partial_forder(cset[1:100,], measure="cont")
  p2 <- partial_forder(cset[-(1:100),], measure="cont")
  expect_equal(combine_forder(list(p1, p2)), forder(cset, measure="cont"))
  p1 <- partial_forder(cset[1:100,], measure="erl")
  p2 <- partial_forder(cset[-(1:100),], measure="erl")
  expect_equal(combine_forder(list(p1, p2)), forder(cset, measure="erl"))
})
