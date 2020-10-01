# Return the part of labels(terms(formula.full)) that is not included in formula.reduced
# Doesn't care if formula.reduced contains something not in formula.full
labels_interesting <- function(formula.full, formula.reduced) {
  mf <- attr(terms(formula.full), "factors")
  mr <- attr(terms(formula.reduced), "factors")
  labs <- character()
  if(attr(terms(formula.full), "intercept") > attr(terms(formula.reduced), "intercept"))
    labs <- c("(Intercept)", labs)
  
  if(length(mf) > 0) {
    varorder <- match(row.names(mr), row.names(mf))
    
    newvars <- setdiff(1:nrow(mf), varorder)
    # What do terms >= 2 mean?
    interesting <- sapply(1:ncol(mf), function(i) {
      term <- mf[,i,drop=FALSE]
      if(any(term[newvars] != 0)) return(TRUE)
      term <- term[varorder]
      term[is.na(term)] <- 0 # Variables only in formula.reduced
      for(j in 1:ncol(mr)) {
        if(all((term!=0) == (mr[,j]!=0))) return(FALSE)
      }
      return(TRUE)
    })
    labs <- c(labs, labels(terms(formula.full))[interesting])
  }
  labs
}

# Check that formula.reduced is nested within formula.full and includes something extra.
check_isnested <- function(formula.full, formula.reduced) {
  # Check that the reduced model is nested within the full model
  if(!all(labels(terms(formula.reduced)) %in% labels(terms(formula.full)))) stop("The reduced model includes some extra variables, not in the full model.")
  if(attr(terms(formula.full), "intercept") < attr(terms(formula.reduced), "intercept")) stop("The reduced model includes intercept, but the full model does not.")
  # Check that the full model includes something in addition to the reduced model
  if(all(labels(terms(formula.full)) %in% labels(terms(formula.reduced)))) stop("The full model should not be equal to the reduced model.")
}
