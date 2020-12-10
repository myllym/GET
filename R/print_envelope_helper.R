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
printhelper_method <- function(x, istest, adj) {
  if(inherits(x, c("global_envelope2d", "combined_global_envelope2d")))
    xdim <- " (2d)"
  else
    xdim <- " (1d)"
  if(inherits(x, c("combined_global_envelope", "combined_global_envelope2d")))
    nstep <- paste0(" (", ifelse(attr(x, "nstep")==1, "one", "two"), "-step)")
  else
    nstep <- ""
  if(is.null(attr(x, "alpha")))
    level <- NULL
  else
    level <- 100*(1-attr(x, "alpha"))
  cat(attr(x, "method"), nstep, xdim, ":\n", sep="")
  if(!istest) { # The case of a central region & fboxplot
    if(!is.null(level)) {
      cat(" * ", level, "% central region based on the measure \"", attr(x, "type"), "\"\n", sep="")
      if(!is.null(attr(x, "M"))) {
        o <- order(attr(x, "M"))[1:5]
        cat(" * The 5 most extreme functions based on \"", attr(x, "type"), "\": ",
            paste0(o, " "), "\n", sep="")
      }
    }
  }
  else { # The case of a global envelope test
    cat(" * Based on the measure: \"", attr(x, "type"), "\"\n",
        " * Significance level of the global test: ", attr(x, "alpha"), "\n", sep="")
    if(adj)
      cat(" * Adjusted level of the global test: ", attr(x, "alpha_star"), "\n", sep="")
    if(!is.null(level)) cat(" * ", level, "% global envelope\n",
        " * p-value of the global test: ", attr(x, "p"), sep="")
    if(!is.null(attr(x, "ties")))
      cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else
      cat("\n")
    if(!is.null(attr(x, "p_interval")))
      cat(" * p-interval       : (", attr(x, "p_interval")[1], ", ",
          attr(x, "p_interval")[2],")\n", sep="")
  }
}

printhelper_method_fboxplotextra <- function(x, combined=FALSE) {
  if(!combined) {
    if(!is.null(attr(x, "outliers"))) outnames <- colnames(attr(x, "outliers"))
    else outnames <- "none"
  }
  else {
    if(!is.null(attr(x, "outliers")[[1]])) outnames <- colnames(attr(x, "outliers")[[1]])
    else outnames <- "none"
  }
  cat(" * Expansion factor: ", attr(x, "factor"), "\n",
      " * Outliers: ", paste0(outnames, " "), "\n", sep="")
}

printhelper_contains <- function(x, istest) {
  if(istest) n <- "global envelope"
  else n <- "central region "
  if('r' %in% names(x)) {
    cat("$r - Argument values                       : ")
    str(x[['r']])
  }
  if('x' %in% names(x)) {
    cat("$x - x-coordinates of pixels               : ")
    str(x[['x']])
    cat("$y - y-coordinates of pixels               : ")
    str(x[['y']])
    cat("$width - Width of pixels at (x, y)         : ")
    str(x[['width']])
    cat("$height - Height of pixels at (x, y)       : ")
    str(x[['height']])
  }
  if('xmin' %in% names(x)) {
    cat("$xmin - Corner coordinates of the pixels   : ")
    str(x[['xmin']])
    cat("$xmax - Corner coordinates of the pixels   : ")
    str(x[['xmax']])
    cat("$ymin - Corner coordinates of the pixels   : ")
    str(x[['ymin']])
    cat("$ymax - Corner coordinates of the pixels   : ")
    str(x[['ymax']])
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
    cat("$lo - Lower boundary of the ", n, ": ", sep="")
    str(x[['lo']])
  }
  if('hi' %in% names(x)) {
    cat("$hi - Upper boundary of the ", n, ": ", sep="")
    str(x[['hi']])
  }
  if('whisker.lo' %in% names(x)) {
    cat("$whisker.lo - Lower boundary of the boxplot: ")
    str(x[['whisker.lo']])
    cat("$whisher.hi - Upper boundary of the boxplot: ")
    str(x[['whisker.hi']])
  }
}

printhelper_contains_combined <- function(x) {
  cat("The object contains a list of ", length(x), " components\n",
      " * each containing: ", paste0("$", names(x[[1]]), " "), "\n", sep="")
}

# A helper function for printing global envelopes
printhelper_ge_base <- function(x, adj=FALSE) {
  if(is.null(attr(x, "p"))) istest <- FALSE
  else istest <- TRUE
  printhelper_method(x, istest, adj)
  if(istest)
    cat(" * Number of r with observed function outside the envelope: ", how_many_outside(x), "\n",
        " * Total number of argument values r                      : ", length(x[['obs']]), "\n", sep="")
  cat("The object contains: \n")
  printhelper_contains(x, istest)
}

# A helper function for printing combined global envelopes
printhelper_ge_combined <- function(x, adj=!is.null(attr(x, "alpha_star"))) {
  if(is.null(attr(x, "p"))) istest <- FALSE
  else istest <- TRUE
  printhelper_method(x, istest, adj)
  printhelper_contains_combined(x)
  if(istest)
    cat(" * Number of r-values with observed function outside the envelope: ",
        paste0(sapply(x, function(x_part) how_many_outside(x_part), simplify=TRUE), " "), "\n",
        " * Total number of argument values r                             : ",
        paste0(sapply(x, function(x_part) length(x_part[['obs']]), simplify=TRUE), " "), "\n", sep="")
}
