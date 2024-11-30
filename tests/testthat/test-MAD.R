test_that("correspondence_MAD_deviation_vs_envelope", {
  cset <- create_curve_set(list(obs=runif(10), sim_m=matrix(runif(10*99), ncol=99)))
  # two-sided unscaled
  res1 <- global_envelope_test(cset, type="unscaled")
  res2 <- deviation_test(cset, scaling = 'none', measure = 'max', savedevs=TRUE)
  expect_equal(attr(res1, "p"), res2$p)
  expect_equal(attr(res1, "M"), res2$devs)
  # one-sided unscaled
  res1 <- global_envelope_test(cset, type="unscaled", alternative="greater")
  res2 <- deviation_test(cset, scaling = 'none', measure = 'max', alternative="greater", savedevs=TRUE)
  expect_equal(attr(res1, "p"), res2$p)
  expect_equal(attr(res1, "M"), res2$devs)
  # st and qdir
  res1 <- global_envelope_test(cset, type="st")
  res2 <- deviation_test(cset, scaling = 'st', measure = 'max', savedevs=TRUE)
  expect_equal(attr(res1, "p"), res2$p)
  expect_equal(attr(res1, "M"), res2$devs)
  res1 <- global_envelope_test(cset, type="qdir")
  res2 <- deviation_test(cset, scaling = 'qdir', measure = 'max', savedevs=TRUE)
  expect_equal(attr(res1, "p"), res2$p)
  expect_equal(attr(res1, "M"), res2$devs)
})
