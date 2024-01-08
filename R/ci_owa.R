ci_owa <- function(x, id, indic_col,atleastjp)
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
  #########                 OWA                 #########
  #######################################################
  ### transpose
  rownames(ci_data) <- x[,id] 
  indic_t <- as.data.frame(t(ci_data))
  
  rows<-n_indic
  cols<-n_unit
  
  ### order each column separately without consider the indicator 
  data_ord <- data.frame(1:n_indic,
                         matrix(NA,ncol=n_unit, nrow=n_indic))
  for(i in 1:n_unit)
  {
    data_ord[,i+1]<- sort(indic_t[,i], decreasing = TRUE)
  }
  colnames(data_ord)<- c(id,colnames(indic_t))
  
  data_ord2<- data_ord[,-1]
  
  #######################################################
  
  n=rows
  alfap=atleastjp/n
  wp <- c(rep(0,n-atleastjp),rep(1/atleastjp,atleastjp))
  wn <- c(rep(1/atleastjp,atleastjp),rep(0,n-atleastjp))
  CI_OWA_p<- matrix(NA, ncol=1, nrow=cols)
  CI_OWA_n<- matrix(NA, ncol=1, nrow=cols)
  
  for(c in 1:cols)
  {
    CI_OWA_n[c,] <- sum(data_ord2[,c]*wp) #more than j
    CI_OWA_p[c,] <- sum(data_ord2[,c]*wn) #at least j
  }
  
  r<- list(CI_OWA_n=CI_OWA_n, CI_OWA_p=CI_OWA_p, wp=wp, wn=wn, ci_method = "owa")
  r$call <- match.call()
  class(r) <- "CI"
  r
}  

