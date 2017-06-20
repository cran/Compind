ci_bod_constr <- function(x,indic_col,up_w,low_w)
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
  

  #library(lpSolve)
  x_num = as.matrix(x_num)

  
  soluz <- rep(0, times=n_unit)
  eff   <- cbind(matrix(0, nrow=n_unit, ncol=(n_indic)))
  final <- cbind(matrix(0, nrow=n_unit, ncol=(n_indic+1)))
  
  for (i in 1:n_unit) 
  {  
    f.obj   <- x_num[i,]
    f.dir   <- c(rep("<=",times=n_unit), rep(">=",times=(n_indic+n_indic)),rep(">=",times=(n_indic)))
    f.rhs   <- c(rep(1,times=n_unit), rep(0,times=(n_indic+n_indic)), rep(0,times=(n_indic)))
    nconstr <- x_num
    Pweight <- cbind(diag(1, nrow=n_indic,ncol=n_indic))
    
    ####### Contraints 
    I = diag(x = 1, n_indic)
    upper = - I+up_w
    lower =  I-low_w
    d3_up  = t(as.vector(f.obj)*upper)
    d3_low = t(as.vector(f.obj)*lower)
    d3 = rbind(d3_up,d3_low)
    ####### End Contraints 
    
    f.con <- rbind(nconstr,d3, Pweight)
    jj <- lp("max", f.obj, f.con, f.dir, f.rhs)
    soluz[i] <- jj$objval
    eff[i,] <- rbind(jj$solution)
    final[i,] <- c(soluz[i], eff[i,])
  }    
  
  w <- ((eff * x_num)/soluz)

  r<-list(ci_bod_constr_est=soluz,
          ci_bod_constr_weights = w, 
          ci_method="bod_constrained")
  r$call<-match.call()
  class(r)<-"CI"
  r
    
  ##return(ci_bod_constr_est)
}


