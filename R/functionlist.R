#' A method function that calculates t-test for each coefficient
#'
#' @param object an object of LReg
#' @param alpha level of significance (default is 0.05)
#' @param ... redundant argument
#'
#' @return estimates, standard error, and confidence intervals of coefficients
#' @export
#'


summary.tconf<-function(object, alpha = 0.05, ...){

  val = object
  for( i in 1:length(val)) assign(names(val)[i], val[[i]])

  #### Calculate SE (t distribution)
  se_t<-
    sqrt(tcrossprod(diag(solve(crossprod(x,x))), sigma2_cor))|>
    `colnames<-`("se")|>
    `rownames<-`(colnames(x))

  lower = (beta_hat - stats::qt(1-alpha/2, n-p-1)*se_t)|>
    `colnames<-`("lower")
  upper = (beta_hat + stats::qt(1-alpha/2, n-p-1)*se_t)|>
    `colnames<-`("upper")

  tidy = cbind(beta_hat, se_t, lower, upper)|>
    as.data.frame()|>
    tibble::rownames_to_column("term")|>
    dplyr::as_tibble()

  (tidy<- dplyr::mutate(tidy, significance = ifelse(0 <=lower | 0 >upper, "[+]","[-]")))

  obj_name<- deparse(substitute(object)) # get name of object and transfer it to string
  obj_env<- pryr::where(obj_name) # get the environment of the object
  object$tconf<- tidy # revise the class of object
  assign(obj_name, object, envir = obj_env) # assign it to object in its original environment

  message('[+] indicates significant, and [-] is insignificant')
  return(tidy)
}

#########################################################################################

#' Calculate the 95% confidence interval for the mean response and future response
#'
#' @param object an object of LReg
#' @param x_star a matrix of new observations
#' @param alpha level of significance (default is 0.05)
#' @param ... redundant argument
#'
#' @return CI of mean response and future response
#' @export
#'

summary.CIresponse = function(object, x_star = NULL, alpha = 0.05, ...){

  if(is.null(x_star)) x_star = object$x[,-1, drop = F]
  if(!is.matrix(x_star)) stop("x_star must be a matrix")
  if(ncol(x_star) != object$p) stop(paste("the column number must equal", object$p))


  x =  cbind(1, x_star) # add intercept

  y_star = x%*%object$beta_hat

  se_mean_res = sqrt(diag(x%*%object$xtx_inv%*%t(x))%*%object$sigma2_cor)
  mean_res<-
    mapply(function(x,y) x + c(-1,1)*y, y_star, stats::qt(1-alpha/2, object$n - object$p-1)*se_mean_res)|>
    t()|>
    `colnames<-`(c("lower_response", "upper_response"))


  se_mean_new_res = sqrt(diag(x%*%object$xtx_inv%*%t(x)+1)%*%object$sigma2_cor)
  mean_new_res<-
    mapply(function(x,y) x + c(-1,1)*y, y_star, stats::qt(1-alpha/2, object$n - object$p-1)*se_mean_new_res)|>
    t()|>
    `colnames<-`(c("lower_future", "upper_future"))

  (result<- cbind(y_star = y_star, mean_res, mean_new_res)|> dplyr::as_tibble())

  obj_name<- deparse(substitute(object)) # get name of object and transfer it to string
  obj_env<- pryr::where(obj_name) # get the environment of the object
  object$CIresponse<- result # revise the class of object
  assign(obj_name, object, envir = obj_env) # assign it to object in its original environment

  return(result)
}

#########################################################################################

#' Calculate joint hypothesis testing for the model parameters
#'
#' @param object an object of LReg
#' @param beta the values of testing
#' @param alpha level of significance (default is 0.05)
#' @param ... redundant argument
#'
#' @return F statistic, lower and upper bound of F distribution at the level of alpha, and significance status
#' @export
#'

summary.JCRpars = function(object, beta = rep(0,object$p+1), alpha = 0.05,...){


  d = object$beta_hat - as.matrix(beta)

  Fstat = (t(d)%*%solve(object$xtx_inv)%*%d)/(object$p+1)/object$sigma2_cor
  ci = stats::qf(c(alpha/2,1-alpha/2), object$p+1, object$n - object$p-1)

  (result<-
  dplyr::tibble(Fstat = c(Fstat), "lower" = ci[1], "upper" = ci[2])|>
    dplyr::mutate(significance = dplyr::if_else(Fstat<= lower | Fstat > upper, "[+]", "[-]")))

  obj_name<- deparse(substitute(object)) # get name of object and transfer it to string
  obj_env<- pryr::where(obj_name) # get the environment of the object
  object$JCRpars<- result # revise the class of object
  assign(obj_name, object, envir = obj_env) # assign it to object in its original environment

  message(paste('[+] indicates significant, and [-] is insignificant.\n alpha =', alpha))
  return(result)
}

#########################################################################################


#' Calculate F statistic for contrast comparisons
#'
#' @param object an object of LReg
#' @param R a matrix of contrast
#' @param r the values of testing (default is 0)
#' @param alpha level of significance (default is 0.05)
#'
#' @return F statistic, lower and upper bound of F distribution at the level of alpha, and significance status
#' @export
#'

ContrastTest = function(object, R, r = rep(0,nrow(R)), alpha = 0.05){

  if(is.vector(R)) stop("R must be in a metrix form")
  if(ncol(R) != (object$p+1)) stop(paste("The number columns of R must equal", object$p+1))
  if(length(r) != nrow(R)) stop(paste("The length of r must equal", nrow(R)))

  r = as.matrix(r)
  k = nrow(r)

  d = R%*%object$beta_hat - r
  v = solve(R%*%object$xtx_inv%*%t(R)*c(object$sigma2_cor))

  Fstat = t(d)%*%v%*%d/k
  ci = stats::qf(c(alpha/2,1-alpha/2), k, object$n - object$p-1)
  (result<-
  dplyr::tibble('Fstat' = c(Fstat), 'lower' = ci[1], 'upper' = ci[2])|>
    dplyr::mutate(significance = dplyr::if_else(Fstat<= lower | Fstat > upper, "[+]", "[-]")))

  message(paste('[+] indicates significant, and [-] is insignificant.\n alpha =', alpha))
  return(result)
}

#########################################################################################
utils::globalVariables(c("x", "sigma2_cor", "beta_hat", "n", "p", "lower", "upper"))





