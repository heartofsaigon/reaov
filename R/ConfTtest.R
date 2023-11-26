#' A function that calculates t-test for each coefficient
#'
#' @param object an object of LReg
#'
#' @return estimates, standard error, and confidence intervals of coefficients
#'
#' @export
#'
#' @examples ConfTtest(LReg(Species ~ Area +Elevation+ Nearest+ Scruz+ Adjacent, mydata1))


ConfTtest = function(object){

  val = object

  for( i in 1:length(val)) assign(names(val)[i], val[[i]])


  #### Calculate SE (t distribution)
  se_t<-
    sqrt(tcrossprod(diag(solve(crossprod(x,x))), sigma2_cor))|>
    `colnames<-`("se")|>
    `rownames<-`(colnames(x))

  lower = (beta_hat - stats::qt(0.975, n-p-1)*se_t)|>
    `colnames<-`("lower")
  upper = (beta_hat + stats::qt(0.975, n-p-1)*se_t)|>
    `colnames<-`("upper")

  tidy = cbind(beta_hat, se_t, lower, upper)|>
    as.data.frame()|>
    tibble::rownames_to_column("term")|>
    dplyr::as_tibble()

  tidy<- dplyr::mutate(tidy, significance = ifelse(0 <=lower | 0 >upper, "yes","no"))
  print(tidy)
  return(append(list(tidy = tidy), val))
}

utils::globalVariables(c("x", "sigma2_cor", "beta_hat", "n", "p"))
