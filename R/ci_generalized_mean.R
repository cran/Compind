ci_generalized_mean <- function(x, indic_col, p, na.rm=TRUE)
{
 # library(psych)
  x_num   = x[,indic_col]
  n_indic <- dim(x_num)[2]
  n_unit  <- dim(x_num)[1]
  
  # Numeric check
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(x_num[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    }
  }   

  
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
  
    
  
  t_x = t(x_num)
  
  if (!na.rm) {
    x_elev_p = x_num^p
    sum_p = apply(x_elev_p,1,sum)
    ci_generalized_mean_est = (1/n_indic*sum_p)^1/p
    
  } else {
    ci_generalized_mean_est <- 1 #### CAMBIARE
  }
  
  r<-list(ci_generalized_mean_est=ci_generalized_mean_est, ci_method="generalized_mean")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)
  
}

