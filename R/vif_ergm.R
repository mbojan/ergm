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



# Unique model term name within a model
# trm = element of ergm_model$terms
.term_name <- function(trm) {
  if(length(trm$coef.names) == 1) return(trm$coef.names)
  l <- strsplit(trm$coef.names, "\\.")
  paste(Reduce(intersect, l), collapse=".")
}

vif_ergm_new <- function(object, drop_terms = NULL, drop_coef = "edges") {
  stopifnot(inherits(object, "ergm"))
  
  em <- ergm_model(object$formula, nw = object$network)
  if(length(em$terms) < 2) stop("model contains fewer than 2 terms")
  
  # DF of terms with names and coef.names
  term_db <- data.frame(
    name = vapply(em$terms, "[[", character(1), "name"),
    term_name = vapply(em$terms, .term_name, character(1)),
    coef.names = lapply(em$terms, "[[", "coef.names"),
    stringsAsFactors = FALSE
  )

  v <- vcov(object)

  # Coef names to be dropped
  if(!is.null(drop_coef)) {
    bad_coef_name <- !(drop_coef %in% em$coef.names)
    if(any(bad_coef_name))
      stop("coefs not in the model: ", 
           paste(drop_coef[bad_coef_name], collapse=", "))
  }
  # Get coef names of terms to be dropped
  if(!is.null(drop_terms)) {
    i.term <- match(drop_terms, term_db$term_name)
    bad_term_name <- is.na(i.term)
    if(any(bad_term_name)) 
      stop("terms not in the model: ", 
           paste(drop_terms[bad_term_name], collapse=", "), "\n", 
           "Available term names: ", paste(term_db$term_name, collapse=", "))
    drop_coef <- c(drop_coef, unlist(term_db$coef.names[i.term]))
  }
  if(!is.null(drop_coef)) {
    i <- match(drop_coef, colnames(v))
    v <- v[-i, -i]
  }
  
  R <- cov2cor(v)
  detR <- det(R)
  
  # Per coefficient
  vifs <- numeric(ncol(R))
  vifs_df <- rep(1, ncol(R))
  for(i in seq(along=colnames(R))) {
    gvec <- R[-i, i, drop=FALSE]
    xcor <- solve(R[-i, -i, drop=FALSE])
    vifs[i] <- 1 / (1 - t(gvec) %*% xcor %*% gvec)
  }
  
  structure(
    list(
      vif = vifs,
      vif_df = vifs_df
    ),
    class = "ergm_vif"
  )
}
