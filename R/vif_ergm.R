#' Variance Inflation Factors for ERGMs
#' 
#' @param object model object as returned by [ergm()]
#' 
#' @return 
#' 
#' @references 
#' Duxbury, S. W. (2018). Diagnosing multicollinearity in exponential random
#' graph models. Sociological Methods & Research <doi:10.1177/0049124118782543>
#' 
#' @export


# Unique model term labels
# 
# object = ergm_model
.term_labels <- function(object) {
  stopifnot(inherits(object, "ergm_model"))
  # Just deparse the term calls as used in the model formula
  vapply(object$terms, function(x) deparse(x$call), character(1))
}



vif_ergm_new <- function(object, drop_terms = NULL, drop_coef = NULL, ...) {
  stopifnot(inherits(object, "ergm"))
  verbose <- getOption("verbose", FALSE)
  
  # Create `ergm_model` object to query the terms
  em <- ergm_model(object$formula, nw = object$network)

  # DF of terms with names, labels, and coef.names
  coefs_per_term <- vapply(em$terms, function(k) length(k$coef.names), integer(1))
  term_db <- data.frame(
    name = rep(vapply(em$terms, "[[", character(1), "name"), coefs_per_term),
    term_label =rep(.term_labels(em), coefs_per_term),
    coef_name = unlist(lapply(em$terms, "[[", "coef.names")),
    stringsAsFactors = FALSE
  )

  # Which rows/cols of the vcov matrix should be dropped?
  # Coefs:
  if(!is.null(drop_coef)) {
    bad_coef_name <- !(drop_coef %in% term_db$coef_name)
    if(any(bad_coef_name))
      stop("coefs not in the model: ", 
           paste(drop_coef[bad_coef_name], collapse=", "))
  }
  # Term labels (expecting deparsed formula summands):
  if(!is.null(drop_terms)) {
    bad_term_label <- !(drop_terms %in% term_db$term_label)
    if(any(bad_term_label)) 
      stop("terms not in the model: ", 
           paste(drop_terms[bad_term_label], collapse=", "), "\n", 
           "Terms in the model: ", paste(term_db$term_label, collapse=", "))
  }
  # We are dropping union of `drop_coefs` and coefs corresponding to
  # `drop_terms`
  drop_coef <- with(
    term_db,
    coef_name[(coef_name %in% drop_coef) | (term_label %in% drop_terms) ]
  )

  # VCov matrix
  v <- vcov(object, ...)
  
  if(length(drop_coef) > 0) {
    i <- match(drop_coef, colnames(v))
    v <- v[-i, -i]
    term_db <- term_db[-i,]
  }
  term_map <- with(term_db, match(term_label, unique(term_label)))
  
  R <- cov2cor(v)
  detR <- det(R)
  
  # Calculate VIFs based on the correlation matrix and term matching vector
  # R = correlation matrix between coefficients
  # a = term assignment matching terms to row/cols of R
  .do_vif <- function(R, a) {
    stopifnot(identical(length(a), nrow(R)))
    stopifnot(identical(length(a), ncol(R)))
    detR <- det(R)
    vif <- numeric(max(a))
    df <- integer(max(a))
    for(term in seq(1, max(a))) {
      i <- which(a == term)
      vif[term] <- det(R[i, i, drop=FALSE]) * det(R[-i, -i, drop=FALSE]) / detR
      df[term] <- length(i)
    }
    data.frame(
      vif = vif, 
      df = df,
      vif_sqrt = vif^(1/(2 * df))
    )
  }
  

  # Return
  structure(
    list(
      vif_coef = cbind(
        coef = term_db$coef_name, 
        .do_vif(R, seq(1, nrow(R)))
      ),
      vif_term = cbind(
        term = unique(term_db$term_label),
        .do_vif(R, term_map)
      )
    ),
    class = "ergm_vif"
  )
}


