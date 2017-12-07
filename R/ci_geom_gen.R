ci_geom_gen <- function(x,indic_col,meth,up_w,low_w,bench)
{

  x_num   = x[,indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  rowProds <- function(X){ apply(X,1,FUN="prod") }
  
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
  
  
  
  if (meth=="EQUAL")  {
    
    geom_BoD = rowProds(x_num ^ (1/n_indic) )
    
    r<-list(ci_mean_geom_est=geom_BoD, ci_method="mean_geom")
    r$call<-match.call()
    class(r)<-"CI"
    r
    
  }
  
  
  if (meth=="BOD")  {  
    
  app = ci_bod_constr(x_num,indic_col,up_w,low_w)
  w = app$ci_bod_constr_weights

  
  # Normalising by benchmark unit vector
  x_num_norm = sweep(x_num, 2, t(x_num[bench,]), "/",check.margin=FALSE)
  
  # calculation data_norm * weights
  mult_BoD = rowProds(x_num_norm ^ w)

  r<-list(ci_geom_bod_est=mult_BoD,
          ci_geom_bod_weights = w, 
          ci_method="geometric_bod")
  r$call<-match.call()
  class(r)<-"CI"
  r
}

r
}


