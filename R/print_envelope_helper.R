how_many_outside <- function(x) {
  switch(get_alternative(x),
         two.sided = {
           sum(x[['obs']]<x[['lo']] | x[['obs']]>x[['hi']])
         },
         less = {
           sum(x[['obs']]<x[['lo']])
         },
         greater = {
           sum(x[['obs']]>x[['hi']])
         })
}

# Helper function for printing object with attributes "alpha", "type", "method"
# and optionally "p", "p_interval", "ties", "alpha_star"
printhelper_method <- function(x, istest, ...) {
  if(!istest) { # The case of a central region
    if(inherits(x, c("fboxplot", "combined_fboxplot"))) {
      cat(attr(x, "method"), ":\n",
          " * ", 100*(1-attr(x, "alpha")), "% central region based on the measure \"", attr(x, "type"), "\"\n",
          " * Expansion factor: ", attr(x, "factor"), "\n",
          " * Outliers: ",
          ifelse(is.null(attr(x, "outliers")), "none", colnames(attr(x, "outliers"))), "\n", sep="")
    }
    else {
      o <- order(attr(x, "k"))[1:5]
      cat(100*(1-attr(x, "alpha")), "% central region based on the measure \"", attr(x, "type"), "\"\n", sep="")
      cat(" * The 5 most extreme functions:", o[1], o[2], o[3], o[4], o[5], "\n")
    }
  }
  else { # The case of a global envelope test
    cat(attr(x, "method"), " based on the measure \"", attr(x, "type"), "\"\n",
        " * Significance level of the test: ", attr(x, "alpha"), "\n",
        " * ", 100*(1-attr(x, "alpha")), "% global envelope\n",
        " * p-value of the test: ", attr(x, "p"), sep="")
    if(!is.null(attr(x, "ties")))
      cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else
      cat("\n")
    if(!is.null(attr(x, "p_interval")))
      cat(" * p-interval       : (", attr(x, "p_interval")[1], ", ",
          attr(x, "p_interval")[2],")\n", sep="")
  }
}

printhelper_contains <- function(x, istest) {
  if(istest) n <- "global envelope"
  else n <- "central region"
  if('r' %in% names(x)) {
    cat("$r - Argument values                       : ")
    str(x[['r']])
  }
  if('obs' %in% names(x)) {
    cat("$obs - Observed function                   : ")
    str(x[['obs']])
  }
  if('central' %in% names(x)) {
    cat("$central - Central function                : ")
    str(x[['central']])
  }
  if('lo' %in% names(x)) {
    cat("$lo - Lower boundary of the", n, ": ")
    str(x[['lo']])
  }
  if('hi' %in% names(x)) {
    cat("$hi - Upper boundary of the", n, ": ")
    str(x[['hi']])
  }
  if('whisker.lo' %in% names(x)) {
    cat("$whisker.lo - Lower boundary of the boxplot: ")
    str(x[['whisker.lo']])
    cat("$whisher.hi - Upper boundary of the boxplot: ")
    str(x[['whisker.hi']])
  }
}
