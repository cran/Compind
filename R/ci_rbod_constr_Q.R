ci_rbod_constr_Q <- function (x, indic_col, low_w=0, pref=NULL, M, B, Q=NULL, Q_ord=NULL, bandwidth) 
{
  x_num = x[, indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  ngood=n_indic - 1
  nbad=1
  
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
  # if(n_indic!=(ngood+nbad)){
  #   stop(paste("The number of simple indicators is not equal to the sum of the
  #              number of good and bad outputs, please insert the right number
  #              of simple indicators in x or choose the right number of columns
  #              in indic_col or indicate the right number of good and bad outputs
  #              in ngood and nbad."))
  # }
  # if(ngood==0){
  #   stop(paste("The number of good outputs (ngood) has to be greater than 0."))
  # }
  # if(nbad==0){
  #   stop(paste("The number of bad outputs (nbad) has to be greater than 0."))
  # }
  
  if(!is.null(Q) & is.null(Q_ord)){ 
  Q = as.matrix(Q)
#   number_exogenous <- dim(Q)[2]
#   nq_cont <- dim(Q)[2]
#   nq_ord  <- 0
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    Q_ord = as.matrix(Q_ord)
#     number_exogenous <- dim(Q_ord)[2]
#     nq_cont <- 0
#     nq_ord  <- dim(Q_ord)[2]
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    Q = as.matrix(Q)
    Q_ord = as.matrix(Q_ord)
#     number_exogenous <- dim(Q)[2] + dim(Q_ord)[2] 
#     nq_cont <- dim(Q)[2]
#     nq_ord  <- dim(Q_ord)[2]
  }
  

  x_num = as.matrix(x_num)
  
  g <- ngood
  b <- nbad
  p <- length(pref)

  bw_cx <- bandwidth
  
  Avg <- apply(x_num, 2, mean)
  
  
  score <- rep(0, times = n_unit)
  eff <- rep(0, times = n_unit)
  w <- cbind(matrix(0, nrow = n_unit, ncol = (n_indic)))
  final <- cbind(matrix(0, nrow = n_unit, ncol = (1+n_indic)))
  kerzi <- matrix(nrow = n_unit, ncol=1)
  
  f.sign <- c(rep(-1, g), rep(-1, b)) ################################
  
  for (i in 1:n_unit) {
    
    if(!is.null(pref))
    {
      # display the number of the NUTS II regions under assessment
      #print(i)
      # Collect the data in a dataframe
      if(!is.null(Q) & is.null(Q_ord)){ 
      dat <- data.frame(Q)
      }
      
      if(is.null(Q) & !is.null(Q_ord)){ 
      dat <- data.frame(apply(Q_ord, 2, ordered))
      }
      
      if(!is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q,apply(Q_ord, 2, ordered))
      }
      
      tdata <- dat[i,]
      # Estimate the densities
      kerz <- npksum(bws=t(bw_cx[i,]),txdat=dat, exdat=tdata, return.kernel.weights=TRUE,cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine") 
      
      #kerz <- npudens(bws=t(bw_cx[i,]),cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine",tdat=tdata,edat=dat)
      kerz<-kerz$kw
      kerzi[i,]<-sum(kerz)
      # collect the performance data
      y_aggr <- data.frame(matrix(x_num[,(g+1):(g+b)],ncol=b),matrix(x_num[,1:g],ncol=g),kerz=kerz)
      yk <- y_aggr[i,c(1:(b+g))]
      q_y <- (b+g)
      subsety<-y_aggr
      score_B <- matrix(nrow=B,ncol=1)
      w_B <- cbind(matrix(0, nrow = B, ncol = (n_indic)))
      final_B <- cbind(matrix(0, nrow = B, ncol = (1+n_indic)))
      fl <- seq(1:n_unit)
      k <- 1
      
      for(k in 1:B)
      {
        draw <- sample(fl, size=M, replace=TRUE, prob=subsety$kerz)
        x_num_sample <- x_num[draw,]
        
        # Different cases of ordVWR
        # Case 1: the preference vector has only 2 elements => only ordVWR_type1  
        if(length(pref)==2)
        {
          f.obj <- c(f.sign * x_num[i, ])
          
          f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                     rep(">=",times=(n_indic-1)),
                     rep(">=", times = g+b))
          
          f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                     rep(0,times=(n_indic-1)),
                     rep(0, times = g+b))  
          
          Constr1 <- c(x_num[i,])
          Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
          
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
          M_imp2[,which(colnames(M_imp)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
          ordVWR_type2 <- (I2*Avg)+ M_imp2
          
          '%!in%' <- function(x,y)!('%in%'(x,y))
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp) %!in% pref),]
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp)==pref[1]),]
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp)==pref[length(pref)]),]
          
          
          Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
          Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
          
          f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1)
          
          jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
          score_B[k] <- jj$objval
          w_B[k, ] <- rbind(jj$solution)
          final_B[k, ] <- c(score_B[k], w_B[k, ])
        }
      
      # Case 2: the preference vector has all indicators         
      
      if(length(pref)==n_indic){    
        
        f.obj <- c(f.sign * x_num[i, ])
        
        f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                   rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)),
                   rep(">=", times = g+b))
        
        f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                   rep(0,times=(n_indic-1)), rep(0,times=(p-2)),
                   rep(0, times = g+b))  
        
        Constr1 <- c(x_num[i,])
        Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
        
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
        score_B[k] <- jj$objval
        w_B[k, ] <- rbind(jj$solution)
        final_B[k, ] <- c(score_B[k], w_B[k, ])
      }
      
      # Case 3: the preference vector has a number of elements >2 and < of the number of indicators 
      
      if(length(pref)>2 & length(pref)<n_indic){
        
        f.obj <- c(f.sign * x_num[i, ])
        
        f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                   rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)),
                   rep(">=", times = g+b))
        
        f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                   rep(0,times=(n_indic-1)), rep(0,times=(p-2)),
                   rep(0, times = g+b))  
        
        Constr1 <- c(x_num[i,])
        Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
        
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
        score_B[k] <- jj$objval
        w_B[k, ] <- rbind(jj$solution)
        final_B[k, ] <- c(score_B[k], w_B[k, ])
      }
    }
  }else{
    # display the number of the NUTS II regions under assessment
    #print(i)
    # Collect the data in a dataframe
      if(!is.null(Q) & is.null(Q_ord)){ 
        dat <- data.frame(Q)
      }
      
      if(is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(apply(Q_ord, 2, ordered))
      }
      
      if(!is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q,apply(Q_ord, 2, ordered))
      }
      
    tdata <- dat[i,]
    # Estimate the densities
    kerz <- npksum(bws=t(bw_cx[i,]),txdat=dat, exdat=tdata, return.kernel.weights=TRUE,cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine") 
    
    #kerz <- npudens(bws=t(bw_cx[i,]),cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine",tdat=tdata,edat=dat)
    kerz<-kerz$kw
    kerzi[i,]<-sum(kerz)
    # collect the performance data
    y_aggr <- data.frame(matrix(x_num[,(g+1):(g+b)],ncol=b),matrix(x_num[,1:g],ncol=g),kerz=kerz)
    yk <- y_aggr[i,c(1:(b+g))]
    q_y <- (b+g)
    subsety<-y_aggr
    score_B <- matrix(nrow=B,ncol=1)
    w_B <- cbind(matrix(0, nrow = B, ncol = (n_indic)))
    final_B <- cbind(matrix(0, nrow = B, ncol = (1+n_indic)))
    fl <- seq(1:n_unit)
    k <- 1
    
    for(k in 1:B)
    {
      draw <- sample(fl, size=M, replace=TRUE, prob=subsety$kerz)
      x_num_sample <- x_num[draw,]
      
      f.obj <- c(f.sign * x_num[i, ])
      
      f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                 rep(">=", times = g+b))
      f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                 rep(0, times = g+b))  
      
      Constr1 <- c(x_num[i,])
      Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
      
      I = diag(x = 1, n_indic)
      VWRg = t((I - low_w)*as.vector(Avg)) 
      Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
      Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
      
      f.con <- rbind(Constr1, Constr2, VWRg, Pweight_g, Pweight_b)
      jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
      score_B[k] <- jj$objval
      w_B[k, ] <- rbind(jj$solution)
      final_B[k, ] <- c(score_B[k], w_B[k, ])
    }
  }
   score[i]<- apply(score_B,2,mean) 
   w[i,]<- apply(w_B,2,mean)  
   final[i,]<- apply(final_B,2,mean)   
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
  
  r <- list(ci_rbod_constr_Q_est = eff, ci_rbod_constr_Q_weights = w, 
            ci_rbod_constr_Q_target = target, ci_method = "rbod_constr_Q")
  r$call <- match.call()
  class(r) <- "CI"
  r
}