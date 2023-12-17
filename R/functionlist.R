#########################################################################################
### Methods
#########################################################################################

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
    `colnames<-`(c("lower_mean_response", "upper_mean_response"))


  se_mean_new_res = sqrt(diag(x%*%object$xtx_inv%*%t(x)+1)%*%object$sigma2_cor)
  mean_new_res<-
    mapply(function(x,y) x + c(-1,1)*y, y_star, stats::qt(1-alpha/2, object$n - object$p-1)*se_mean_new_res)|>
    t()|>
    `colnames<-`(c("lower_response", "upper_response"))

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
#' @param alpha level of significance (default is 0.05)
#' @param ... redundant argument
#'
#' @return F statistic, quantile of F distribution at the level of 1-alpha, and significance status
#' @export
#'

summary.JCRpars = function(object, alpha = 0.05, ...){

  beta = rep(0,object$p+1)
  d = object$beta_hat - as.matrix(beta)

  Fstat = (t(d)%*%solve(object$xtx_inv)%*%d)/(object$p+1)/object$sigma2_cor # inverse of inverse is X'X
  ci = stats::qf(c(1-alpha), object$p+1, object$n - object$p-1)

  upper_point = paste0((1-alpha)*100,"%")

  (result<-
  dplyr::tibble(Fstat = c(Fstat), !!rlang::sym(upper_point) := ci)|>
    dplyr::mutate(significance = dplyr::if_else(Fstat > !!rlang::sym(upper_point), "[+]", "[-]")))

  obj_name<- deparse(substitute(object)) # get name of object and transfer it to string
  obj_env<- pryr::where(obj_name) # get the environment of the object
  object$JCRpars<- result # revise the class of object
  assign(obj_name, object, envir = obj_env) # assign it to object in its original environment

  message(paste('[+] indicates significant, and [-] is insignificant.\n alpha =', alpha))
  return(result)
}

#########################################################################################
### functions
#########################################################################################


#' Calculate F statistic for contrast comparisons
#'
#' @param object an object of LReg
#' @param R a matrix of contrast
#' @param alpha level of significance (default is 0.05)
#'
#' @return F statistic, quantile of F distribution at the level of 1-alpha, and significance status
#' @export
#'

ContrastTest = function(object, R, alpha = 0.05){

  if(is.vector(R)) stop("R must be in a matrix form")
  if(ncol(R) != (object$p+1)) stop(paste("The number columns of R must equal", object$p+1))

  r = rep(0,nrow(R))
  r = as.matrix(r)
  k = nrow(r)

  d = R%*%object$beta_hat - r
  v = solve(R%*%object$xtx_inv%*%t(R)*c(object$sigma2_cor)) # bring the sigma_cor to the numerator

  Fstat = t(d)%*%v%*%d/k
  ci = stats::qf(c(1-alpha), k, object$n - object$p-1)

  upper_point = paste0((1-alpha)*100,"%")

  (result<-
  dplyr::tibble('Fstat' = c(Fstat),  !!rlang::sym(upper_point) := ci)|>
    dplyr::mutate(significance = dplyr::if_else(Fstat > !!rlang::sym(upper_point), "[+]", "[-]")))

  message(paste('[+] indicates significant, and [-] is insignificant.\n alpha =', alpha))
  return(result)
}

#########################################################################################

#' Calculate F statistic to test sub-vector of coefficients
#'
#' @param object an object of LReg
#' @param reduced a vector of numbers relating to the ordering of reduced predictors (The ordering depends on the formula in the object of LReg)
#' @param alpha level of significance (default is 0.05)
#'
#' @return F statistic,  quantile of F distribution at the level of 1-alpha, and significance status
#' @export
#'


SubsetBetaTest = function(object, reduced = 1 , alpha = 0.05){

  if(sum(reduced < 1)>0 | sum(reduced> (object$p+1))>0){
    stop(paste("reduced predictors are within", 1, "and", object$p+1))
  }

  fullmodX = crossprod(object$y_hat, object$y_hat)

  X = object$x[, -reduced]
  H0modX2 = c(t(object$y)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%object$y)

  modX1 = fullmodX - H0modX2

  Ftest = c(modX1/length(reduced)/crossprod(object$e_hat, object$e_hat)*(object$n - object$p-1))
  ci = stats::qf(c(1-alpha), length(reduced), (object$n - object$p-1))

  upper_point = paste0((1-alpha)*100,"%")


  (result<-
      dplyr::tibble(`Fstat` = Ftest, !!rlang::sym(upper_point) := ci )|>
      dplyr::mutate(significance = dplyr::if_else(Fstat > !!rlang::sym(upper_point), "[+]", "[-]")))

  message(paste('[+] indicates significant, and [-] is insignificant.\n alpha =', alpha))

  return(result)
}


#########################################################################################
utils::globalVariables(c("x", "sigma2_cor", "beta_hat", "n", "p", "lower", "upper", "Fstat",
                         ":=", "0%"))





