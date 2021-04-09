ci_ampi = function (x, indic_col, gp, time, polarity, penalty = "NEG") 
{
  x_num = x[, indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  if (n_indic < 2) {
    stop(paste("There must be at least two simple indicators!"))
  }
  for (i in seq(1, n_indic)) {
    if (!is.numeric(x_num[, i])) {
      stop(paste("Data set not numeric at column:", 
                 i))
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
  timef = levels(as.factor(time))
  ci_ampi_est <- matrix(0, nrow = (n_unit/length(timef)), ncol = length(timef))
  for (t in 1:length(timef)) {
    x_num_t <- x_num[time == timef[t], ]
    ci_norm <- x_num_t
    max <- matrix(0, nrow = 1, ncol = n_indic)
    min <- matrix(0, nrow = 1, ncol = n_indic)
    delta <- matrix(0, nrow = 1, ncol = n_indic)
    minim <- matrix(0, nrow = 1, ncol = length(gp))
    maxim <- matrix(0, nrow = 1, ncol = length(gp))
    min = apply(x_num, 2, function(x) min(x, na.rm = TRUE))
    max = apply(x_num, 2, function(x) max(x, na.rm = TRUE))
    delta = (max - min)/2
    maxim = gp + delta
    minim = gp - delta
    
    for (i in seq(1, n_indic)) {
      if (polarity[i] == "POS") {
        ci_norm[, i] = (((x_num_t[, i]) - minim[i])/(maxim[i] - 
                                                       minim[i])) * 60 + 70
      }
      if (polarity[i] == "NEG") {
        ci_norm[, i] = ((maxim[i] - (x_num_t[, i]))/(maxim[i] - 
                                                       minim[i])) * 60 + 70
      }
      if (polarity[i] != "NEG" & polarity[i] != "POS") {
        stop("Please check polarity!")
      }
    }
    Ma_z <- apply(ci_norm, 1, mean)
    Sqm_z <- (apply(ci_norm, 1, function(x) {var(x)*(length(x)-1)/length(x)}))^0.5
    cv = Sqm_z/Ma_z
    ifelse(penalty=="POS", ci_ampi_est[, t] <- round(Ma_z + (Sqm_z * cv), 3), 
           ci_ampi_est[, t] <- round(Ma_z - (Sqm_z * cv), 3))
  }
  ci_ampi_est <- as.data.frame(ci_ampi_est)
  colnames(ci_ampi_est) <- timef
  r <- list(ci_ampi_est = ci_ampi_est, ci_method = "ampi")
  r$call <- match.call()
  class(r) <- "CI"
  return(r)
}
