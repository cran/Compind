plot_M_robust_ci <- function(x,indic_col,method,mvector,B,dir,interval=NULL)
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
  
  uni <- as.matrix(seq(1, 1, len = dim(x_num)[1]))
  ci_data = cbind(x_num, uni)
  
  msuper <- matrix(0, nrow = length(mvector), ncol = 1)
  
  for (m in 1:length(mvector)) {
    if(method=="rbod")
    {
      scores <- ci_rbod(x, indic_col, M=mvector[m],B=B)$ci_rbod_est
      msuper[m,1] <- length(scores[scores >= 1])/nrow(ci_data)
    }
    if(method=="rbod_dir")
    {
      scores <- ci_rbod_dir(x, indic_col, M=mvector[m],B=B, dir=dir)$ci_rbod_dir_est
      msuper[m,1] <- length(scores[scores >= 1])/nrow(ci_data)
    }
    if(method=="rbod_mdir")
    {
      if (is.null(interval)) interval <- 0.05
      scores <- ci_rbod_mdir(x, indic_col, M=mvector[m],B=B,interval=interval)$ci_rbod_mdir_est
      # MDir_RBoD scores are bounded by 1: share of units on the robust frontier
      msuper[m,1] <- length(scores[scores >= 1 - 1e-6])/nrow(ci_data)
    }
  }
  
  if (method == "rbod_mdir") {
    main_lab <- "Percentage of units on the robust frontier by subset size M"
    y_lab    <- "% of units with score = 1"
  } else {
    main_lab <- "Percentage of outperforming units by subset size M"
    y_lab    <- "% of outperforming units"
  }
  plot(x = mvector, y = msuper[, 1], type = "b", lwd = 2, main = main_lab, 
       xlab = c("M"), ylab = y_lab, 
       ylim = c(min(msuper), max(msuper)))
  graphics::lines(x = mvector, y = msuper[, 1], type = "b", 
                  lwd = 2)
  
}