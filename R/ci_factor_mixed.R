

ci_factor_mixed <- function (x, indic_col, method = "ONE", dim = 3) {
  
  #library(FactoMineR)
  x_num   = x[,indic_col]
  n_indic <- length(indic_col)
  n_unit <- nrow(x)

  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (is.na(x_num[i,j]))
      {
        message(paste("Pay attention: NA values at column:",i,", row",j,". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
        #       options(warn=-2)  
      }
    }
  }    
  
  
  
  if (method == "ONE") {
    famd_results <- FAMD(x[, indic_col], ncp = 1, graph = FALSE)
    pesi_fatt = as.matrix(var(famd_results$ind$coord)^2)
    ci_factor_est = famd_results$ind$coord
    r <- list(ci_factor_est = ci_factor_est, loadings_fact = pesi_fatt, 
              ci_method = "factor_mixed")
    
  } else if (method == "ALL") {
    # Make sure that ncp is not greater than min(n_unit-1, n_indic)
    ncp_to_use <- min(n_unit - 1, n_indic)
    famd_results <- FAMD(x[, indic_col], ncp = ncp_to_use, graph = FALSE)
    # Calculate the weights for each factor
    pesi_fatt = apply(famd_results$var$coord^2, 2, sum) / ncp_to_use
    # Compute the composite indicator scores
    ci_factor_est = famd_results$ind$coord %*% pesi_fatt
    r <- list(ci_factor_est = ci_factor_est, loadings_fact = pesi_fatt, 
              ci_method = "factor_mixed")
    
  } else if (method == "CH") {
    famd_results <- FAMD(x[, indic_col], ncp = dim, graph = FALSE)
    pesi_fatt = apply(famd_results$var$coord^2, 2, sum) / n_indic
    pesi_fatt = pesi_fatt[1:dim]
    ci_factor_est = famd_results$ind$coord[, 1:dim] %*% pesi_fatt
    r <- list(ci_factor_est = ci_factor_est, loadings_fact = pesi_fatt, 
              ci_method = "factor_mixed")
  } else {
    stop("Please check method!")
  }
  
  r$call <- match.call()
  class(r) <- "CI"
  return(r)
}


