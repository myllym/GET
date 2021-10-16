how_many_outside <- function(x) {
  switch(get_alternative(x),
         two.sided = {
           sum(x[['obs']]<x[[env_loname(attr(x, "alpha"), all=FALSE)]] | x[['obs']]>x[[env_hiname(attr(x, "alpha"), all=FALSE)]])
         },
         less = {
           sum(x[['obs']]<x[[env_loname(attr(x, "alpha"), all=FALSE)]])
         },
         greater = {
           sum(x[['obs']]>x[[env_hiname(attr(x, "alpha"), all=FALSE)]])
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
      for(i in 1:length(level)) cat(" * ", level[i], "% central region based on the measure \"", attr(x, "type"), "\"\n", sep="")
      if(!is.null(attr(x, "M"))) {
        o <- order(attr(x, "M"))[1:5]
        cat(" * The 5 most extreme functions based on \"", attr(x, "type"), "\": ",
            paste0(o, " "), "\n", sep="")
      }
    }
  }
  else { # The case of a global envelope test
    cat(" * Based on the measure: \"", attr(x, "type"), "\"\n", sep="")
    if(!is.null(level))
      cat(" * ", paste0(paste(level), collapse=", "), "% global envelope", sep="")
    if(adj)
      cat("(adjusted levels ", paste0(paste(100*(1-attr(x, "alpha_star"))), collapse=", "), ")\n", sep="")
    else
      cat("\n")
    cat(" * p-value of the global test: ", attr(x, "p", exact=TRUE), sep="")
    if(!is.null(attr(x, "ties")))
      cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else
      cat("\n")
    if(!is.null(attr(x, "p_interval")))
      cat(" * p-interval       : (", attr(x, "p_interval")[1], ", ",
          attr(x, "p_interval")[2],")\n", sep="")
    cat(" * Significance level of the global test: ", attr(x, "alpha")[1], "\n", sep="")
    if(adj)
      cat(" * Adjusted level of the global test: ", attr(x, "alpha_star")[1], "\n", sep="")
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
  lo <- grep('lo', names(x))
  if(length(lo) == 1) {
    cat("$lo - Lower boundary of the ", n, ": ", sep="")
    str(x[['lo']])
  }
  if(length(lo) > 1) {
    cat(paste0("$", names(x)[lo], collapse=", "), " - Lower boundaries\n", sep="")
  }
  hi <- grep('hi', names(x))
  if(length(hi) == 1) {
    cat("$hi - Upper boundary of the ", n, ": ", sep="")
    str(x[['hi']])
  }
  if(length(hi) > 1) {
    cat(paste0("$", names(x)[hi], collapse=", "), " - Upper boundaries\n", sep="")
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
printhelper_ge_base <- function(x, adj=!is.null(attr(x, "alpha_star"))) {
  if(is.null(attr(x, "p", exact=TRUE))) istest <- FALSE
  else istest <- TRUE
  printhelper_method(x, istest, adj)
  if(istest) {
    cat(" * Number of r with observed function outside the envelope: ", how_many_outside(x), "\n",
        " * Total number of argument values r                      : ", length(x[['obs']]), "\n", sep="")
  }
  cat("The object contains: \n")
  printhelper_contains(x, istest)
}

# A helper function for printing combined global envelopes
printhelper_ge_combined <- function(x, adj=!is.null(attr(x, "alpha_star"))) {
  if(is.null(attr(x, "p", exact=TRUE))) istest <- FALSE
  else istest <- TRUE
  printhelper_method(x, istest, adj)
  printhelper_contains_combined(x)
  if(istest)
    cat(" * Number of r-values with observed function outside the envelope: ",
        paste0(sapply(x, function(x_part) how_many_outside(x_part), simplify=TRUE), " "), "\n",
        " * Total number of argument values r                             : ",
        paste0(sapply(x, function(x_part) length(x_part[['obs']]), simplify=TRUE), " "), "\n", sep="")
}
