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

vif_ergm <- function(my.ergm){
  require(ergm)
  
  cor.mat<-cov2cor(my.ergm$covar) #calculate correlation matrix
  corr5<-cor.mat[-c(1),-c(1)] ##omit edges term
  
  corr5<-corr5[!is.na(corr5[1:nrow(corr5)]),]
  corr5<-corr5[,which(!is.na(corr5[1,1:ncol(corr5)]))]
  
  VIFS<-matrix(0,nr=1,nc=ncol(corr5))
  
  for(i in 1:ncol(corr5)){
    
    gvec<-as.vector(corr5[-c(i),i]) ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec<-t(gvec)            
    xcor<-solve(corr5[-c(i),-c(i)]) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq<-tgvec%*%xcor%*%gvec
    VIFS[1,i]<-1/(1-Rsq)
  }
  ##Columns are covariates as they appear in the SAOM object
  colnames(VIFS)<-names(my.ergm$coef[-c(1)])
  message("Higher values indicate greater correlation.\nVIF > 20 is concerning, VIF > 100 indicates severe collinearity.")
  VIFS
}