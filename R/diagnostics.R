#########################################################################################
#diagnostic plots
#########################################################################################

#' Return a list of diagnostic plots
#' @param object an object of LReg
#' @return list of diagnostic plots
#' @export
#'

Diagnostic_Plots = function(object){
  diagnostics_tib <-  dplyr::as_tibble(c(list(
    e_hat = object$e_hat,
    e_hat_st = object$e_hat/diag(object$hat_mat),
    y_hat = object$y_hat),
    data.frame(object$x[,2:dim(object$x)[2]])
  ))
  x_names <- colnames(data.frame(object$x[,2:dim(object$x)[2]]))

  #make an empty list of plots, names are already filled in
  list_names <- c("residuals_v_fitted",
                  "residuals_hist",
                  paste("residuals_v_", x_names, sep = ''), #residuals_v_xi
                  "qq_standardized_residuals"
  )
  list_plots <- stats::setNames(vector(mode = "list", length=length(list_names)), list_names)

  #residuals vs fitted values
  list_plots$residuals_v_fitted <- diagnostics_tib |>
    ggplot2::ggplot(ggplot2::aes(x = object$y_hat, y = object$e_hat)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "red")+
    ggplot2::labs(title = "Residuals vs Fitted Values",
         x = "Fitted Value",
         y = "Residual")

  #residuals histogram
  list_plots$residuals_hist <- diagnostics_tib |>
    ggplot2::ggplot(ggplot2::aes(x = e_hat)) +
    ggplot2::geom_histogram()+
    ggplot2::labs(title = "Residuals Histogram",
         x = "Residual")

  #residuals vs covariates
  for(x_name in x_names){
    list_plots[[paste("residuals_v_", x_name, sep = '')]] <- diagnostics_tib |>
      ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x_name), y = e_hat)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "red")+
      ggplot2::labs(title = paste("Residuals vs", x_name),
           x = x_name,
           y = "Residual")
  }

  #qq plot of standardized residuals vs theoretical quantile
  list_plots$qq_standardized_residuals <- diagnostics_tib |>
    ggplot2::ggplot(ggplot2::aes(sample = e_hat_st))+
    ggplot2::geom_qq() +
    ggplot2::geom_qq_line(color = "red") +
    ggplot2::labs(title = "Normal QQ Plot Standardized Residuals",
         x = "Theoretical",
         y = "Standardized Residual")

  return(list_plots)
}

utils::globalVariables(c("e_hat", "e_hat_st"))
