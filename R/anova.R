#' A function that calculated the R squared and adjusted R squared values of a
#' model
#'
#' @param object an object of LReg
#'
#' @return R squared and adjusted R squared value of object
#' @export
#'

Rsquared = function(object){

  val1 = 1-sum(object$e_hat**2)/sum((object$y-mean(object$y))**2)

  n = length(object$y)
  p = dim(object$x)[2]
  val2 = 1-(1-val1)*n/(n-p)

  R = list(R_sq = val1, adj_R_sq = val2)

  return(R)
}


#' A method function that calculates the ANOVA table given a vector of responses
#' an a vector of categorical covariate observation
#'
#' @param response vector of responses
#' @param covariates vector of categorical covariate observation
#'
#' @return factor and residual degrees of freedom, between group sum of squares,
#' total sum of squares, mean squared errors, F-statistic, p-value
#' @export
#'
ano2 = function(response, covariates){

  response = unlist(response)
  covariates = unlist(covariates)

  y_bar = mean(response)
  new_x = as.factor(covariates)
  levels(new_x) = 1:length(levels(new_x))

  g = 1
  i = 1
  count = 0

  y_gbar = rep(0,length(levels(new_x)))
  while(g<length(levels(new_x))+1){
    while (i<length(covariates)+1){
      if(new_x[i]==g){
        count = count + 1
        y_gbar[g]=y_gbar[g]+response[i]
      }
      i = i+1
    }
    if(count>0){
      y_gbar[g]=y_gbar[g]/count
    }
    i=1
    g = g+1
    count = 0
  }

  BSS =0
  WSS =0
  g = 1
  i = 1

  while(g<length(levels(new_x))+1){
    while (i<length(covariates)+1){
      if(new_x[i]==g){
        WSS = WSS + (response[i]-y_gbar[g])**2
        BSS = BSS + (y_gbar[g]-y_bar)**2
      }
      i = i+1
    }
    i = 1
    g = g + 1
  }

  n = length(covariates)
  G = length(levels(new_x))
  n_p = n-G

  G_p = G-1
  F_stat = (BSS/(G-1))*(n-G/WSS)
  p_value = stats::pf(F_stat,G_p,n_p, lower.tail=FALSE)

  results = list(factor_DG = G-1, residual_DG = n-G, BSS = BSS,BSS_MSE = BSS/n,
       TSS = WSS+BSS, TSS_MSE = (WSS+BSS)/n, F_statistic = F_stat,
       p_value = p_value)

  return(results)
}




