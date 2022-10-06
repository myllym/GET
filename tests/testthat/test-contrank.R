test_that("contrank", {
  for(n in c(1,2,3,4,5,6,10, 100)) {
    x <- rnorm(n)
    r <- rank(x, ties.method = "average")
    cr <- contrank(x)
    expect_equal(all(r-1 <= cr & cr <= r), TRUE)
    x <- sample(1:3, n, replace=TRUE)
    r <- rank(x, ties.method = "average")
    cr <- contrank(x)
    expect_equal(all(r-1 <= cr & cr <= r), TRUE)
  }
})
