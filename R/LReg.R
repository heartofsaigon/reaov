#' A main function that calculates useful values used to fit the linear model
#'
#' @param myformula a formula of the model that includes response and covariates
#' @param mydata data where we fit the linear model
#' @param myclass assign a class to the object
#'
#' @return a list of useful values
#'
#' @export
#'
#' @examples LReg(Species ~ Area +Elevation+ Nearest+ Scruz+ Adjacent, mydata1)

LReg = function(myformula, mydata, myclass = "tconf"){

  fo = as.character(myformula)
  response = fo[2]
  covariate = fo[3]|>stringr::str_extract_all("[:alpha:]+")|>
    unlist()
  covariate<- c("1", covariate)

  intercept = rep(1, nrow(mydata))
  y = mydata[,response]|> as.matrix();
  x = cbind(mydata,`1` = intercept)[,covariate, drop = FALSE]|>
    as.matrix()

  n = nrow(x); p = ncol(x)-1

  xtx_inv = solve(crossprod(x,x))
  hat_mat = tcrossprod(x%*%xtx_inv,x)
  i_mat = diag(rep(1, n))

  y_hat = hat_mat%*%y
  e_hat = (i_mat - hat_mat)%*%y

  beta_hat = (xtx_inv%*%crossprod(x,y))|>
    `colnames<-`("estimate")

  sigma2_naive = crossprod(e_hat,e_hat)/n
  sigma2_cor = sigma2_naive*n/(n-p-1)

  result<-
  list(
    y = y, x = x,
    n = n, p = p,
    hat_mat = hat_mat, i_mat = i_mat, xtx_inv = xtx_inv,
    y_hat = y_hat, e_hat = e_hat,
    beta_hat = beta_hat,
    sigma2_naive = sigma2_naive, sigma2_cor = sigma2_cor
  )

  class(result)<- c(myclass, class(result))
  result
}

