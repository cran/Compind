ci_geom_bod_intertemp <- function(x0,x1,indic_col,up_w,low_w,bench)
{

  Indic_t0   = x0[,indic_col]
  Indic_t1   = x1[,indic_col]
  n_indic <- dim(as.matrix(Indic_t0))[2]
  n_unit <- dim(as.matrix(Indic_t0))[1]
  lower_w = low_w
  upper_w = up_w  
  
  ################################################################################
  #                                                                              #
  #  Intertemporal setting: period t0-t1                                         #
  #                                                                              #
  ################################################################################
  

  rowProds <- function(X){ apply(X,1,FUN="prod") }
  
  # Normalising by benchmark unit vector
  Indic_norm_t0 = sweep(Indic_t0, 2, t(Indic_t0[bench,]), "/")
  Indic_norm_t1 = sweep(Indic_t1, 2, t(Indic_t1[bench,]), "/")
  
  # Constrained BoD weights
  CI_t0 = ci_bod_constr(Indic_t0, indic_col,upper_w,lower_w)
  w_t0  = CI_t0$ci_bod_constr_weights
  CI_t1 = ci_bod_constr(Indic_t1, indic_col,upper_w,lower_w)
  w_t1  = CI_t1$ci_bod_constr_weights
  
  
  # Index calculation ########################################################################
  
  # Computations of the Indicator Change Effect (Norm) A
  ChangeEffect_Norm_t1_vs_t0 <- rowProds((Indic_norm_t1/Indic_norm_t0)^((w_t1+w_t0)/2))
  
  # Computations of the Indicator Change Effect (NO Norm) B
  ChangeEffect_t1_vs_t0      <- rowProds((Indic_t1/Indic_t0)^((w_t1+w_t0)/2))
  
  # Computations of the Benchmark Change Effect C
  benc_m = matrix(0, nrow=n_unit, ncol=n_indic)
  for (j in 1:n_unit){
    benc_m[j,] <- (Indic_t1[bench,]/Indic_t0[bench,])^(-(w_t1[j,]+w_t0[j,])/2)
  }
  BenchmarkEffect_t1_vs_t0   <- rowProds(benc_m)
  
  # Computations of the Weight Change Effect D
  WeightEffect_t1_vs_t0    <- rowProds((Indic_norm_t1*Indic_norm_t0)^((w_t1-w_t0)/2))
  
  # Computation of the Overall Change
  OverallChange_t1_vs_t0     <- #ChangeEffect_Norm_t1_vs_t0 * BenchmarkEffect_t1_vs_t0 
                                ChangeEffect_t1_vs_t0 * BenchmarkEffect_t1_vs_t0 # WeightEffect_t1_vs_t0
  
  #############################################################################################  
  
  # cbind
  Intertemp_t1_vs_t0 <- cbind(OverallChange_t1_vs_t0, ChangeEffect_t1_vs_t0,
                              BenchmarkEffect_t1_vs_t0, WeightEffect_t1_vs_t0) 
  Intertemp_t1_vs_t0  
  

  r<-list(ci_geom_bod_intertemp_est=Intertemp_t1_vs_t0,
          ci_method="Intertemporal_effects_Geometric_BoD")
  r$call<-match.call()
  class(r)<-"CI"
  r


}


