#' Variance Inflation Factors for ERGMs
#' 
#' @param model model object as returned by [ergm()]
#' 
#' @return 
#' 
#' @author 
#' Based on original code by Scott W. Duxbury
#' 
#' @references 
#' Duxbury, S. W. (2018). Diagnosing multicollinearity in exponential random
#' graph models. Sociological Methods & Research <doi:10.1177/0049124118782543>
#' 
#' @export

vif_ergm <- function(object, ...){
  # Correlation matrix without the edges term
  # TODO this should come from a set of simulated networks
  cormat <- stats::cov2cor(object$covar)[-1, -1]

  # MB: Not sure when there would be NAs in `covar`. Offsets perhaps? Commenting
  # for now.
  # 
  # corr5<-corr5[!is.na(corr5[1:nrow(corr5)]),]
  # corr5<-corr5[,which(!is.na(corr5[1,1:ncol(corr5)]))]
  
  
  VIFS <- numeric(ncol(cormat))
  
  for(i in 1:ncol(cormat)){
    # Correlations of i-th with the rest
    gvec <- cormat[-i , i, drop=FALSE] ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec <- t(gvec)            
    xcor <- solve( cormat[-i, -i] ) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq <- tgvec %*% xcor %*% gvec
    VIFS[i] <- 1/(1-Rsq)
  }
  
  names(VIFS) <- names(object$coef[-1])
  #  message("Higher values indicate greater correlation.\nVIF > 20 is concerning, VIF > 100 indicates severe collinearity.")
  VIFS
}



# Unique model term labels
# 
# object = ergm_model
.term_labels <- function(object) {
  stopifnot(inherits(object, "ergm_model"))
  # Just deparse the term calls as used in the model formula
  vapply(object$terms, function(x) deparse(x$call), character(1))
}



vif_ergm_new <- function(object, drop_terms = NULL, drop_coef = NULL) {
  stopifnot(inherits(object, "ergm"))
  
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

  # Get vcov matrix
  v <- vcov(object)

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
    # We will drop union of drop_coefs and coefs corresponding to drop_terms
    drop_coef <- unique(c(drop_coef, term_db$coef_name[term_db$term_label %in% drop_terms]))
  }
  # Drop selected rowcols from vcov
  if(!is.null(drop_coef)) {
    i <- match(drop_coef, colnames(v))
    v <- v[-i, -i]
  }
  
  R <- cov2cor(v)
  detR <- det(R)
  
  # Per-coefficient VIFs
  vifs <- numeric(ncol(R))
  vifs_df <- rep(1, ncol(R))
  for(i in seq(along=colnames(R))) {
    gvec <- R[-i, i, drop=FALSE]
    xcor <- solve(R[-i, -i, drop=FALSE])
    vifs[i] <- 1 / (1 - t(gvec) %*% xcor %*% gvec)
  }
  
  # TODO: Per-term VIFs
  
  # Return
  structure(
    list(
      vif_coef = data.frame(
        name = rownames(R),
        VIF = vifs,
        df = vifs_df
      )
    ),
    class = "ergm_vif"
  )
}
