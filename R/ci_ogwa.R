ci_ogwa <- function(x, id, indic_col, atleastjp, coords,
                    kernel = "bisquare", adaptive = F, bw, 
                    p = 2, theta = 0, longlat = F, dMat)
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
  ci_data = as.matrix(cbind(x_num))
  
  
  #######################################################
  #########                 OGWA                 ########
  #######################################################
  
  n=n_indic
  alfap=atleastjp/n
  wp <- c(rep(0,n-atleastjp),rep(1/atleastjp,atleastjp))
  wn <- c(rep(1/atleastjp,atleastjp),rep(0,n-atleastjp))
  
  ### matrix of all the distances between all the GW model calibration points and data points
  dist<- matrix(NA,nrow=n_unit,ncol=n_unit)
  
  for(i in 1:n_unit)
  {
    dist[,i]<-gw.dist(dp.locat=coords, focus=i, longlat=longlat) 
  }
  
  #dist<-gw.dist(dp.locat=coords,longlat=TRUE) 
  
  data_sp <- SpatialPointsDataFrame(coords, data.frame(ci_data))
  
  bw <- bw.gwss.average(data_sp,
                        vars=colnames(data_sp@data), 
                        kernel=kernel,
                        adaptive=adaptive, 
                        longlat=longlat,
                        dMat=dist)
  
  CI_OGWA_n <- matrix(NA,nrow=n_unit, ncol=n_unit)
  CI_OGWA_p <- matrix(NA,nrow=n_unit, ncol=n_unit)
  ci_data_gw.vi <- matrix(NA,nrow=n_unit, ncol=n_indic)
  
  for (i in 1:n_unit) {
    dist.vi <- gw.dist(dp.locat=coords, focus=i, longlat=longlat)
    wt.vi <- gw.weight(dist.vi, bw=mean(bw[1,]), kernel=kernel,
                       adaptive=adaptive)
    #use<- wt.vi>0.5
    for(c in 1:n_indic)
    {
      ci_data_gw.vi[,c] <- ci_data[,c]*wt.vi
    }
    
    indic_t_gw.vi <- as.data.frame(t(ci_data_gw.vi))
    
    rows<-n_indic
    cols<-dim(ci_data_gw.vi)[1]
    
    ### order each column separately without consider the indicator 
    data_ord_gw.vi <- data.frame(1:n_indic,
                                 matrix(NA,ncol=dim(ci_data_gw.vi)[1], nrow=n_indic))
    for(j in 1:dim(ci_data_gw.vi)[1])
    {
      data_ord_gw.vi[,j+1]<- sort(indic_t_gw.vi[,j], decreasing = TRUE)
    }
    colnames(data_ord_gw.vi)<- c(id,colnames(indic_t_gw.vi))
    
    data_ord_gw2.vi<- data_ord_gw.vi[,-1]
    
    #######################################################
    
    CI_OGWA_p.vi<- matrix(NA, ncol=1, nrow=cols)
    CI_OGWA_n.vi<- matrix(NA, ncol=1, nrow=cols)
    
    for(c in 1:cols)
    {
      CI_OGWA_n.vi[c,] <- sum(data_ord_gw2.vi[,c]*wp) #more than j
      CI_OGWA_p.vi[c,] <- sum(data_ord_gw2.vi[,c]*wn) #at least j
    }
    
    CI_OGWA_n[,i] <- CI_OGWA_n.vi
    CI_OGWA_p[,i] <- CI_OGWA_p.vi
    
  }
  
  CI_OGWA_n_mean=rowMeans(CI_OGWA_n)
  CI_OGWA_p_mean=rowMeans(CI_OGWA_p)
  
  r<- list(CI_OGWA_n=CI_OGWA_n_mean, CI_OGWA_p=CI_OGWA_p_mean, wp=wp, wn=wn,
           bw=bw,
           ci_method = "ogwa")
  r$call <- match.call()
  class(r) <- "CI"
  r
}  

