#' Variance Inflation Factors for ERGMs
#' 
#' @param object ERGM model specificaton as returned by [ergm_model()] or a
#'   fitted model model object as returned by [ergm()]. See Methods section
#'   below.
#' @param ... Other arguments passed to/from other methods.
#' 
#' @return An object of class `vif_ergm`...
#' 
#' @references Duxbury, S. W. (2018). Diagnosing multicollinearity in
#'   exponential random graph models. *Sociological Methods & Research*
#'   \doi{10.1177/0049124118782543}.
#'
#'   Fox, J., Monette, G. (1992). Generalized Collinearity Diagnostics. *Journal
#'   of the American Statistical Association* 87(417):178--183
#'   \doi{10.2307/2290467}.
#' 
#' @export
vif_ergm <- function(object, ...) UseMethod("vif_ergm")





#' @describeIn vif_ergm Based on `ergm_model` object and VCov matrix
#' 
#' @param drop_terms character vector of term labels
#' @param drop_coef character vector of coefficient names
#' @param v Variance-covariance matrix
#' 
#' @export
vif_ergm.ergm_model <- function(object, v, drop_terms=NULL, drop_coef = NULL) {
  # DF of terms with names, labels, and coef.names
  coefs_per_term <- vapply(object$terms, function(k) length(k$coef.names), integer(1))
  term_db <- data.frame(
    name = rep(vapply(em$terms, "[[", character(1), "name"), coefs_per_term),
    term_label =rep(.term_labels(object), coefs_per_term),
    coef_name = unlist(lapply(object$terms, "[[", "coef.names")),
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
  
  if(length(drop_coef) > 0) {
    i <- match(drop_coef, colnames(v))
    v <- v[-i, -i]
    term_db <- term_db[-i,]
  }
  term_map <- with(term_db, match(term_label, unique(term_label)))
  
  R <- cov2cor(v)
  detR <- det(R)
  
  # Calculate VIFs based on the correlation matrix and term matching vector.
  # Based on Fox & Monette (1992).
  # 
  # R = correlation matrix between coefficients based on vcov(model)
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
    class = "vif_ergm"
  )
}

#' @describeIn vif_ergm Method for [ergm()] modles
#' 
#' @export
vif_ergm.ergm <- function(object, ...) {
  # Create `ergm_model` object to query the terms
  em <- ergm_model(object$formula, nw = object$network)
  vif_ergm.ergm_model(em, ...)
}


# Unique model term labels
# 
# object = ergm_model
.term_labels <- function(object) {
  stopifnot(inherits(object, "ergm_model"))
  # Just deparse the term calls as used in the model formula
  vapply(object$terms, function(x) deparse(x$call), character(1))
}
