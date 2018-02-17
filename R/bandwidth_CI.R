bandwidth_CI <- function(x, indic_col, ngood, nbad, Q=NULL, Q_ord=NULL)
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

  
  x_num = as.matrix(x_num)
  
  g <- ngood
  b <- nbad
  
  if(!is.null(Q) & is.null(Q_ord)){ 
    Q = as.matrix(Q)
    number_exogenous <- dim(Q)[2]
    nq_cont <- dim(Q)[2]
    nq_ord  <- 0
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    Q_ord = as.matrix(Q_ord)
    number_exogenous <- dim(Q_ord)[2]
    nq_cont <- 0
    nq_ord  <- dim(Q_ord)[2]
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    Q = as.matrix(Q)
    Q_ord = as.matrix(Q_ord)
    number_exogenous <- dim(Q)[2] + dim(Q_ord)[2] 
    nq_cont <- dim(Q)[2]
    nq_ord  <- dim(Q_ord)[2]
  }
  
  # Create (n x (b+g)) matrix with only FALSE-values
flag_ii<-matrix(FALSE,nrow=length(x_num[,1]),ncol=(b+g))
# Flag all data values which satisfy the conditions below with TRUE-values

flag_gi <- apply(matrix(x_num[,1:g],ncol=g), 2, function(x) x<=median(x))
flag_bi <- apply(matrix(x_num[,(g+1):(g+b)],ncol=b), 2, function(x) x>=median(x))
flag_ii <- cbind(flag_gi,flag_bi)

# Flag the NUTS II regions which satisfy the (b+g) conditions
flagg <- matrix(nrow=length(x_num[,1]),ncol=1)
for (j in (1:length(x_num[,1])))
{
  flagg[j] <- all(flag_ii[j,])
}
# Flag for the observations who satisfy all criteria the data and exogenous variables with TRUE-values 
if(!is.null(Q) & is.null(Q_ord)){ 
  flag_QQ <- subset(Q,subset=flagg,drop=TRUE)   # the continuous exogenous variables
  data_fr1 <- data.frame(flag_QQ)
}

if(is.null(Q) & !is.null(Q_ord)){ 
  flag_ordd <- subset(Q_ord,subset=flagg,drop=TRUE)  # the discrete (ordered) exogenous variables
  data_fr1 <- data.frame(apply(flag_ordd, 2, ordered))
}

if(!is.null(Q) & !is.null(Q_ord)){ 
  flag_QQ <- subset(Q,subset=flagg,drop=TRUE)   # the continuous exogenous variables
  flag_ordd <- subset(Q_ord,subset=flagg,drop=TRUE)  # the discrete (ordered) exogenous variables
  # Collect the concerned data in a dataframe
  data_fr1 <- data.frame(flag_QQ,apply(flag_ordd, 2, ordered))
}

flag_yy <- subset(x_num,subset=flagg,drop=TRUE)   # the performance variables

# Compute the bandwiths using the performance data and exogenous data
bw1 <-  npcdensbw(ydat=flag_yy,xdat=data_fr1,cykertype="epanechnikov",cxkertype="epanechnikov",bwmethod="cv.ls",oxkertype="liracine",itmax = 1000)   ##set itmax to 1000 nmulti to 4 and tol high
# Take the bandwiths for the exogenous variables
bwfull <- bw1$xbw

# Create an empty (nx1)-matrix
aha<- matrix(nrow=length(x_num[,1]),ncol=1)
# Create an empty matrix to store the endogenous bandwiths for the exogenous variables
bw_cx <- matrix(nrow=length(x_num[,1]),ncol=number_exogenous)
# Set i equal to the start value 1 (start loop with DMU 1)
i <- 1
# Define the loop to indicate which NUTS II regions dominate region i
for (i in (1:length(x_num[,1])))
{
  #Define a (n x (b+g)) matrix with all FALSE values
  flag_i<-matrix(FALSE,nrow=length(x_num[,1]),ncol=(b+g))
  # Flag the data values that dominate the values of DMU i
  flag_g <- apply(matrix(x_num[,1:g],ncol=g), 2, function(x) x<=x[i])
  flag_b <- apply(matrix(x_num[,(g+1):(g+b)],ncol=b), 2, function(x) x>=x[i])
  flag_i <- cbind(flag_g,flag_b)
  #Define a (n x 1) matrix 
  flag <- matrix(nrow=length(x_num[,1]),ncol=1)
  # Define a loop to indicate which NUTS II-regions dominate region i
  for (j in (1:length(x_num[,1])))
  {
    flag[j] <- all(flag_i[j,])
  }
  # print(i)
  # Take the subset of data values that correspond to the flagged NUTS II-regions 

  # Flag for the observations who satisfy all criteria the data and exogenous variables with TRUE-values 
  if(!is.null(Q) & is.null(Q_ord)){ 
    flag_Q <- matrix(subset(Q,subset=flag,drop=TRUE), ncol=nq_cont)   # the continuous exogenous variables
    # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
    data_fr <- data.frame(flag_Q)
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    flag_ord <- matrix(subset(Q_ord,subset=flag,drop=TRUE), ncol=nq_ord)  # the discrete (ordered) exogenous variables
    # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
    data_fr <- data.frame(apply(flag_ord, 2, ordered))
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    flag_Q <- matrix(subset(Q,subset=flag,drop=TRUE), ncol=nq_cont)   # the continuous exogenous variables
    flag_ord <-  matrix(subset(Q_ord,subset=flag,drop=TRUE), ncol=nq_ord)  # the discrete (ordered) exogenous variables
    # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
    data_fr <- data.frame(flag_Q,apply(flag_ord, 2, ordered))
  }
  
  flag_y <- subset(x_num,subset=flag,drop=TRUE)      # the performance variables
  
  #print(sum(flag))     # display the number of NUTS-regions dominating region i
  aha[i]=sum(flag)
  # Define a loop to determine the proper bandwith based on the number of dominating regions
  if((sum(flag))<=40)               # if number of dominating regions lower than or equal to threshold
  {
    bw_cx[i,] <-bwfull              # bandwith equal to what was derived endogenously above
  } else
  {                               # if number of dominating regions higher than threshold

    # Compute the bandwiths using the performance data and exogenous data
    bw <-  npcdensbw(ydat=flag_y,xdat=data_fr,cykertype="epanechnikov",cxkertype="epanechnikov",bwmethod="cv.ls",oxkertype="liracine",itmax = 10000,nmulti=4)   ##set itmax to 1000 nmulti to 4 and tol high
    # Take the bandwiths for the exogenous variables
    bw_cx[i,] <- bw$xbw
    #      print(bw_cx[i,])
  }
}

# save the output of the first step
#write.matrix(bw_cx, file = "welbandwidth.txt", sep ="\t")

# before opening the .txt file, first reset the options for system separators
# Open the file with the endogenously estimated bandwiths (for the exogenous variables)
#bw_cx<-read.table("welbandwidth.txt", sep = "\t")
bw_cx2=bw_cx
# Adjust the estimated bandwiths for the discrete (ordered) exogenous variables
if(!is.null(Q) & is.null(Q_ord)){ 
  bw_cx2 <- bw_cx
}

if(is.null(Q) & !is.null(Q_ord)){ 
  bw_cx2[apply(bw_cx, 2, function(x) x>0.9999999)]<-0.9999
}

if(!is.null(Q) & !is.null(Q_ord)){ 
  bw_cx2[,nq_cont+1:(nq_cont+nq_ord-2)][apply(bw_cx[,nq_cont+1:(nq_cont+nq_ord-2)], 2, function(x) x>0.9999999)]<-0.9999
}

# Specify the final version of the matrix with endogenous bandwiths
bw_cx<-bw_cx2

#write.matrix(bw_cx, file = "welbandwidth_adj.txt", sep ="\t")

r <- list(bandwidth = bw_cx,ci_method = "bandwidth_CI")
r$call <- match.call()
class(r) <- "CI"
r

}
