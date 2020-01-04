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
