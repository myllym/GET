test_that("contrank", {
  expectrank <- function(x) {
    r <- rank(x, ties.method = "average")
    cr <- contrank(x)
    expect_equal(all(r-1 <= cr & cr <= r), TRUE)
  }
  expectrank(c(1,2,2))
  expectrank(c(2,2,3))
  for(n in c(1,2,3,4,5,6,10, 100)) {
    expectrank(rnorm(n))
    expectrank(sample(1:3, n, replace=TRUE))
  }
})
