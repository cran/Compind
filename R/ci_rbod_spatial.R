
ci_rbod_spatial <- function(x, indic_col, M=20, B=100, W) 
{
  
  x_num   = x[,indic_col]
  n_indic = dim(x_num)[2]
  n_unit  = dim(x_num)[1]

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
  
  ###################################
  
  W = as.matrix(W)
  #library(spdep)
  nb <- mat2listw(W, row.names = NULL, style="M")$neighbours
  
  uni <- as.matrix(seq(1, 1, len = n_unit))
  data_per_spRbod = cbind(x_num,uni)
  
  dim = dim(x)[1]
  eff_cond = matrix(NA, nrow=dim, ncol=1)
  
  
  for (i in (1:dim)) {
    
    dati_front_cond = data_per_spRbod[c(nb[[i]]),]
    dati_order_riga = data_per_spRbod[i,]
    rbod_cond = try(orderm(base = dati_order_riga, 
                           frontier = dati_front_cond, 
                           noutput = n_indic, 
                           orientation=2, M = M, B = B))
    if ('try-error' %in% class(rbod_cond)) next
    eff_cond[i] = 1/rbod_cond$eff
  }
  r <- list(ci_rbod_spatial_est = eff_cond, 
            ci_method = "rbod_spatial")
}

