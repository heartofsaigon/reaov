#########################################################################################
#Model Comparison Test
#########################################################################################

#' Return the f statistic for the comparison of a larger vs a smaller model, null hypothesis is small model is correct
#' @param large_model a formula
#' @param small_model a formula, y must be same as large_model and all variables must also be in large_model
#' @return f statistic
#' @export
#'


ModelComparison <- function(large_model, small_model){
  #no error checking, but should be added to see if models are included properly
  del_p <- large_model$p - small_model$p
  p_small <- small_model$p
  n <- large_model$n

  large_RSS <- sum(large_model$e_hat**2)
  small_RSS <-  sum(small_model$e_hat**2)

  F_test <- ((small_RSS - large_RSS)/del_p )/(large_RSS/(n-p_small))
  p_value <- stats::pf(F_test, del_p, n-p_small,lower.tail = FALSE)
  result <- list(
    F_stat = F_test,
    p_value = p_value
  )
  return(result)
}
