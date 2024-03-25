
ci_rbod_mdir <- function(x, indic_col, M, B, interval=0.05)
{

  confidence_interval <- function(vector, interval) {

    vec_sd <- sd(vector, na.rm = TRUE)
    n <- length(vector[!is.na(vector)])
    vec_mean <- mean(vector, na.rm = TRUE)
    error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
    result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
    return(result)
  }
  

  x_num   = x[,indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit  <- dim(as.matrix(x_num))[1]
  
  # Numeric check
  if (n_unit<=M)
  {
    stop("M is greater than or equal to the number of total units")
  }

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

  ########################################################    
  ## Costruzione matrici
  ci_mdir_est       <- matrix(nrow=B,ncol=M)
  ci_mdir_est_exact <- matrix(nrow=B,ncol=n_unit)
  conf <- matrix(nrow=2,ncol=n_unit)
  CImeaSPEC_tot     <- matrix(0,nrow=n_unit,ncol=n_indic)
  CIdirecttot       <- matrix(0,nrow=n_unit,ncol=n_indic)
  list_dir1 <- c()
  list_dir2 <- c()  
  
  
  ########################################################    
  # ciclo per B
  for(k in 1:B)
  {
    
    ci_data = as.matrix(x_num[sample(nrow(x_num), M), ])
    units = as.numeric(row.names(ci_data))
    uni <- as.matrix(seq(1, 1, len = dim(ci_data)[1]))
    den_CImeaSPEC = cbind(matrix(0, nrow = M, ncol = (n_indic)))
    CImeaSPEC = cbind(matrix(0, nrow = M, ncol = (n_indic)))
    CI <- mea(uni,ci_data,  ORIENTATION = "out", RTS="crs")
    
    for (j in seq(1,n_indic)) 
    {
      den_CImeaSPEC[,j] <- ci_data[,j]+((CI$eff)*CI$direct[,j])
      CImeaSPEC[,j]     <- ci_data[,j]/den_CImeaSPEC[,j]
    }
    
    ci_mdir_est[k,] <- 1-rowSums((CI$direct*CI$eff)/rowSums(ci_data+((CI$eff)*CI$direct)))
    
    
    ########################################################    
    # metto i risultati nelle giuste posizioni
    for(z in 1:M)
    {
      ci_mdir_est_exact[k,units[z]] = ci_mdir_est[k,z]
      for (j in seq(1,n_indic)) 
      {
        CImeaSPEC_tot[units[z],j] = CImeaSPEC[z,j]
        CIdirecttot[units[z],j]   = CI$direct[z,j]
      }
      
      list_dir1[[k]] <- CImeaSPEC_tot
      list_dir2[[k]] <- CIdirecttot
      
    }
    cat(paste('\014', 'Loop:',k, "of total", B))
  }
  
  
  ########################################################    
  ###### Calcolo misure finali
  
  E_ci_mdir_est = colMeans(ci_mdir_est_exact, na.rm = TRUE)
  for(w in 1:n_unit) {
    conf[,w] = confidence_interval(ci_mdir_est_exact[,w], 0.95)
  }
  CImeaSPEC_tot = mat_mean(list_dir1, na.rm = TRUE)
  CIdirecttot   = mat_mean(list_dir2, na.rm = TRUE)


  ########################################################    
  # output list

  r<-list(ci_rbod_mdir_est  = E_ci_mdir_est,
          conf  = t(rbind(lower_ci = conf[1,],
                          upper_ci = conf[2,])),
          ci_rbod_mdir_spec = CImeaSPEC_tot, 
          ci_rbod_mdir_dir  = CIdirecttot, 
          ci_method="rbod_mdir")
  r$call<-match.call()
  class(r)<-"CI"
  r
}


