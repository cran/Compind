
ci_smaa_constr <- function(x,indic_col, rep, label, low_w=NULL)
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
  
  # Constraint check
  if (!is.null(low_w)) {
  n_low_w <- dim(as.matrix(low_w))[1]
      if (n_indic != n_low_w)
    {
      stop("Lower bound vector must be set as a vector of the same size as the number of simple indicators")
      }
  }


    ############# Creazione dati                ################
  ############# meas	= Criteria measurements ################
  
  N <- rep        # numero di ripetizioni 
  m <- n_unit     # alternatives = numero soggetti da valutare
  n <- n_indic    # criteria     = numero di indicatori semplici
  
  arr = array(0, dim=c(N,m,n))
  
  for(k in 1:N){
    for(j in 1:n){
      arr[k,,j] = x_num[,j]
    }
  }
  

  ######## Creazione pesi [che sommano a 1] ###################
  ######## pref	= Weights                   ###################
  
  pesi = array(0, dim=c(N,n))
  
  if (is.null(low_w)) {
    for(j in 1:n){
      pesi[,j] = runif(N, min = 0, max=1)
    }
  }

  if (!is.null(low_w)) {
    for(j in 1:n){
      pesi[,j] = runif(N, min = low_w[j], max=1)
    }
  }

  if (!is.null(low_w)) {
    somma_low = sum(low_w)
    # Numeric check
    if (somma_low>1)
    {
      stop(paste("Lower constraint too high. Lower constraint * Number of indicators=",somma_low))
    }  
  }
    
  # Normalizzo la somma dei w a 1 
  for(j in 1:n){
    pesi[,j] =pesi[,j]/sum(pesi[,j])
  }
  


  
  ############################################################
  ########### SMAA ###########################################
  
  # Costruzione della rank frequency table
  values <- smaa.values(meas=arr, pref=pesi)
  ranks  <- smaa.ranks(values)
  
  ranks_freq= array(0, dim=c(m,m))
  for(i in 1:N){
    for(j in 1:m){
      ranks_freq[ranks[i,j],j]  = ranks_freq[ranks[i,j],j] + 1
    }
  }
  
  colnames(ranks_freq) = label
  ranks_freq_t = t(ranks_freq)

  ################# rango medio ###########################
  
  somma <- matrix(0, nrow = m, ncol = 1)
  for(i in 1:m){
    for(j in 1:m){
      somma[i] = somma[i] + (j * ranks_freq_t[i,j]/rep)
    }             
  }
  rownames(somma) = label
  colnames(somma) = "Average_rank"
  


r<-list(ci_smaa_constr_rank_freq=ranks_freq_t,
        ci_smaa_constr_values = values,
        ci_smaa_constr_average_rank = somma,
        ci_method="ci_smaa_constr")
r$call<-match.call()
class(r)<-"CI"
r

}


