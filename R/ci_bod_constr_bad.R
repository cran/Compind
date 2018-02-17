ci_bod_constr_bad <- function(x, indic_col, ngood=1, nbad=1, low_w=0, pref=NULL) 
{
  x_num = x[, indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  for (i in seq(1, n_indic)) {
    if (!is.numeric(x_num[, i])) {
      stop(paste("Data set not numeric at column:", i))
    }
  }
  for (i in seq(1, n_unit)) {
    for (j in seq(1, n_indic)) {
      if (is.na(x_num[i, j])) {
        message(paste("Pay attention: NA values at column:", 
                      i, ", row", j, ". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
      }
    }
  }
  if(!is.null(pref)){
      if (!is.character(pref)) {
       stop(paste("The preference vector (pref) is not of character type"))
      }
  }
  if(n_indic!=(ngood+nbad)){
    stop(paste("The number of simple indicators is not equal to the sum of the
               number of good and bad outputs, please insert the right number
               of simple indicators in x or choose the right number of columns
               in indic_col or indicate the right number of good and bad outputs
               in ngood and nbad."))
    }
  if(ngood==0){
    stop(paste("The number of good outputs (ngood) has to be greater than 0."))
  }
  if(nbad==0){
    stop(paste("The number of bad outputs (nbad) has to be greater than 0."))
  }
  
  
  x_num = as.matrix(x_num)
  
  g <- ngood
  b <- nbad
  p <- length(pref)

  Avg <- apply(x_num, 2, mean)
  
  score <- rep(0, times = n_unit)
  eff <- rep(0, times = n_unit)
  w <- cbind(matrix(0, nrow = n_unit, ncol = (n_indic)))
  final <- cbind(matrix(0, nrow = n_unit, ncol = (1+n_indic)))
  
  f.sign <- c(rep(-1, g), rep(1, b))

  for (i in 1:n_unit) {
    
    if(!is.null(pref))
    {
      
  # Different cases of ordVWR
    # Case 1: the preference vector has only 2 elements => only ordVWR_type1  
      if(length(pref)==2)
      {
             f.obj <- c(f.sign * x_num[i, ])
             
             f.dir <- c("==",rep(">=", times = n_unit), rep(">=", times = g+b), 
                        rep(">=",times=(n_indic-1)),
                        rep(">=", times = g+b))
             
             f.rhs <- c(1, rep(0, times = n_unit), rep(0, times = g+b), 
                        rep(0,times=(n_indic-1)),
                        rep(0, times = g+b))  
             
             Constr1 <- c(x_num[i,])
             Constr2 <- cbind(t(f.sign * t(x_num)))
             
             I = diag(x = 1, n_indic)
             VWRg = t((I - low_w)*as.vector(Avg)) 
             
             ### OrdVWR_type1
             I = diag(x = 1, n_indic)
             colnames(I) <- colnames(x_num)
             I[,pref[1]] <- 0
             
             M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
             colnames(M_imp) <- colnames(x_num)
             M_imp[,pref[1]] <- 1
             
             ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
             ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
             
             Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
             Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
             
             f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1)
             jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
             score[i] <- jj$objval
             w[i, ] <- rbind(jj$solution)
             final[i, ] <- c(score[i], w[i, ])
             
     # Case 2: the preference vector has all indicators         
      }
      
      if(length(pref)==n_indic){    
       
            f.obj <- c(f.sign * x_num[i, ])
            
            f.dir <- c("==",rep(">=", times = n_unit), rep(">=", times = g+b), 
                       rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)), 
                       rep(">=", times = g+b))
            
            f.rhs <- c(1, rep(0, times = n_unit), rep(0, times = g+b), 
                       rep(0,times=(n_indic-1)), rep(0, times = (p-2)),
                       rep(0, times = g+b))  
            
            Constr1 <- c(x_num[i,])
            Constr2 <- cbind(t(f.sign * t(x_num))) 
            
            I = diag(x = 1, n_indic)
            VWRg = t((I - low_w)*as.vector(Avg)) 
            
            ### OrdVWR_type1
            I = diag(x = 1, n_indic)
            colnames(I) <- colnames(x_num)
            I[,pref[1]] <- 0
            
            M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
            colnames(M_imp) <- colnames(x_num)
            M_imp[,pref[1]] <- 1
            
            ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
            ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
            
            ### OrdVWR_type2
            I2 = diag(x = 1, n_indic)
            colnames(I2) <- colnames(x_num)
            I2[,pref[length(pref)]] <- 0
            
            M_imp2 =  matrix(0,nrow=(n_indic), ncol=n_indic)
            colnames(M_imp2) <- colnames(x_num)
            rownames(M_imp2) <- colnames(x_num)
            M_imp2[,which(colnames(M_imp2)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
            ordVWR_type2 <- (I2*Avg)+ M_imp2
            colnames(ordVWR_type2) <- colnames(x_num)
            rownames(ordVWR_type2) <- colnames(x_num)
            
            #'%!in%' <- function(x,y)!('%in%'(x,y))
            #ordVWR_type2 <- ordVWR_type2[-which(colnames(ordVWR_type2) %!in% pref),]
            ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[length(pref)]),]
            ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[1]),]
            
            Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
            Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
            
            f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1, ordVWR_type2)
            jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
            score[i] <- jj$objval
            w[i, ] <- rbind(jj$solution)
            final[i, ] <- c(score[i], w[i, ])
            
    # Case 3: the preference vector has a number of elements >2 and < of the number of indicators 
    }
      if(length(pref)>2 & length(pref)<n_indic){
      
            f.obj <- c(f.sign * x_num[i, ])
            
            f.dir <- c("==",rep(">=", times = n_unit), rep(">=", times = g+b), 
                       rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)), 
                       rep(">=", times = g+b))
            
            f.rhs <- c(1, rep(0, times = n_unit), rep(0, times = g+b), 
                       rep(0,times=(n_indic-1)), rep(0, times = (p-2)),
                       rep(0, times = g+b))  
            
            Constr1 <- c(x_num[i,])
            Constr2 <- cbind(t(f.sign * t(x_num))) 
            
            I = diag(x = 1, n_indic)
            VWRg = t((I - low_w)*as.vector(Avg)) 
            
            ### OrdVWR_type1
            I = diag(x = 1, n_indic)
            colnames(I) <- colnames(x_num)
            I[,pref[1]] <- 0
            
            M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
            colnames(M_imp) <- colnames(x_num)
            M_imp[,pref[1]] <- 1
            
            ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
            ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
            
            ### OrdVWR_type2
            I2 = diag(x = 1, n_indic)
            colnames(I2) <- colnames(x_num)
            I2[,pref[length(pref)]] <- 0
            
            M_imp2 =  matrix(0,nrow=(n_indic), ncol=n_indic)
            colnames(M_imp2) <- colnames(x_num)
            rownames(M_imp2) <- colnames(x_num)
            M_imp2[,which(colnames(M_imp2)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
            ordVWR_type2 <- (I2*Avg)+ M_imp2
            colnames(ordVWR_type2) <- colnames(x_num)
            rownames(ordVWR_type2) <- colnames(x_num)
            
            '%!in%' <- function(x,y)!('%in%'(x,y))
            ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2) %!in% pref),]
            ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[length(pref)]),]
            ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[1]),]
            
            Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
            Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
            
            f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1, ordVWR_type2)
            jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
            score[i] <- jj$objval
            w[i, ] <- rbind(jj$solution)
            final[i, ] <- c(score[i], w[i, ])
    }
 }else{
      f.obj <- c(f.sign * x_num[i, ])
      
      f.dir <- c("==",rep(">=", times = n_unit), rep(">=", times = g+b), 
                 rep(">=", times = g+b))
      f.rhs <- c(1, rep(0, times = n_unit), rep(0, times = g+b), 
                 rep(0, times = g+b))  
      
      Constr1 <- c(x_num[i,])
      Constr2 <- cbind(t(f.sign * t(x_num)))
      
      I = diag(x = 1, n_indic)
      VWRg = t((I - low_w)*as.vector(Avg)) 
      Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
      Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
      
      f.con <- rbind(Constr1, Constr2, VWRg, Pweight_g, Pweight_b)
      jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
      score[i] <- jj$objval
      w[i, ] <- rbind(jj$solution)
      final[i, ] <- c(score[i], w[i, ])
    }
  }
  eff <- (1/(1+score))
  
  colnames(w) <- colnames(x_num)
  
  # Write the composite scores and the indicator target values in one result matrix
  
  target <- matrix(nrow=n_unit,ncol=(b+g)) 
  for(j in 1:g)
  {
  target[,j] <- x_num[,j] + (Avg[j]*score)
  }
  for(j in (g+1):(g+b))
  {
    target[,j] <- x_num[,j] - (Avg[j]*score)
  }
  colnames(target) <- colnames(x_num)
  
  r <- list(ci_bod_constr_bad_est = eff, ci_bod_constr_bad_weights = w, 
            ci_bod_constr_bad_target = target, ci_method = "bod_constr_bad")
  r$call <- match.call()
  class(r) <- "CI"
  r
}