ci_bod_mdir <- function(x,indic_col)
{
  x_num   = x[,indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
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
  
  ci_data = as.matrix(cbind(x_num))
  uni <- as.matrix(seq(1, 1, len = dim(ci_data)[1]))
  den_CImeaSPEC = cbind(matrix(0, nrow = n_unit, ncol = (n_indic)))
  CImeaSPEC = cbind(matrix(0, nrow = n_unit, ncol = (n_indic)))
    
  CI <- mea(uni,ci_data,  ORIENTATION = "out", RTS="crs")
  
  for (j in seq(1,n_indic)) 
  {
  den_CImeaSPEC[,j] <-ci_data[,j]+((CI$eff)*CI$direct[,j])
  CImeaSPEC[,j] <- ci_data[,j]/den_CImeaSPEC[,j]
  }
  
  ci_mdir_est <- 1-rowSums((CI$direct*CI$eff)/rowSums(ci_data+((CI$eff)*CI$direct)))
  
  r<-list(ci_bod_mdir_est=ci_mdir_est,ci_bod_mdir_spec=CImeaSPEC, ci_bod_mdir_dir = CI$direct, ci_method="bod_mdir")
  r$call<-match.call()
  class(r)<-"CI"
  r
}






