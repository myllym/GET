test_that("stepdown", {
  expect_equal(determine_ordering(c(1,1,1,3,1,1,3,3), "erl"), c(2,1,5,3,6,7,4,8))
  expect_equal(determine_ordering(c(3,3,1,3,3,4,3,3,3), "erl"), c(3,2,4,8,1,7,5,9,6))
  expect_equal(determine_ordering(c(0.3, 0.5, 1.2, 0.2, 1.3), "area"), c(1,4,2,3,5))

  expectrankorder <- function(x, m) {
    # Test that ordering agree with rank ordering
    R <- ceiling(x)
    n <- length(x)
    o <- determine_ordering(x, m)
    expect_false(is.unsorted(R[o]))
  }
  for(n in seq(2, 100, length=10)) {
    for(m in seq(1, 20, length=10)) {
      expectrankorder(1+round(runif(n, 0, m)), "erl")
      expectrankorder(runif(n, 0, m), "area")
    }
  }
})
