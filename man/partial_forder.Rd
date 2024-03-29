% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forder.r
\name{partial_forder}
\alias{partial_forder}
\alias{combine_forder}
\title{Functional ordering in parts}
\usage{
partial_forder(
  curve_set,
  measure = c("erl", "rank", "cont", "area"),
  alternative = c("two.sided", "less", "greater")
)

combine_forder(ls)
}
\arguments{
\item{curve_set}{A \code{curve_set} object, usually a part of a larger \code{curve_set}.
(No missing or infinite values allowed.)}

\item{measure}{The measure to use to order the functions from the most extreme to the least extreme
one. Must be one of the following: 'rank', 'erl', 'cont', 'area', 'max', 'int', 'int2'. Default is 'erl'.}

\item{alternative}{A character string specifying the alternative hypothesis.
Must be one of the following: "two.sided" (default), "less" or "greater".
The last two options only available for types \code{'rank'}, \code{'erl'},
\code{'cont'} and \code{'area'}.}

\item{ls}{List of objects returned by partial_forder}
}
\value{
See \code{\link{forder}}
}
\description{
If the functional data doesn't comfortably fit in memory it is possible to
compute functional ordering by splitting the domain of the data (voxels in
a brain image), using \code{partial_forder} on each part and finally
combining the results with \code{combine_forder}.
}
\examples{
data("abide_9002_23")
res <- lapply(list(1:100, 101:200, 201:261), function(part) {
  set.seed(123) # When using partial_forder, all parts must use the same seed.
  fset <- frank.flm(nsim=99, formula.full = Y ~ Group + Sex + Age,
                  formula.reduced = Y ~ Group + Sex,
                  curve_sets = list(Y = abide_9002_23$curve_set[part,]),
                  factors = abide_9002_23$factors, savefuns = "return")
  partial_forder(fset, measure="erl")
})
combine_forder(res)
}
\seealso{
\code{\link{forder}}
}
