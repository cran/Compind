ci_bod_time_decomp <- function(data_long,
                               id_col,
                               time_col,
                               indic_cols) {
  
  stopifnot(is.data.frame(data_long))
  stopifnot(all(c(id_col, time_col, indic_cols) %in% names(data_long)))
  
  df <- data_long
  df <- df[order(df[[time_col]], df[[id_col]]), ]
  
  times <- unique(df[[time_col]])
  times <- sort(times)
  
  if (length(times) < 2) {
    stop("At least two time periods are required to compute the dynamic decomposition.")
  }
  
  if (anyDuplicated(data_long[, c(id_col, time_col)])) {
    stop("Each unit-time combination must appear only once in 'data_long'.")
  }
  
  run_bod <- function(d) {
    X  <- d[, indic_cols, drop = FALSE]
    CI <- ci_bod(X, indic_col = seq_len(ncol(X)))
    as.numeric(CI$ci_bod_est)
  }
  
  # 1) contemporaneous
  E_C <- rep(NA_real_, nrow(df))
  for (tt in times) {
    idx <- which(df[[time_col]] == tt)
    E_C[idx] <- run_bod(df[idx, , drop = FALSE])
  }
  
  # 2) sequential
  E_S <- rep(NA_real_, nrow(df))
  for (tt in times) {
    idx_ref <- which(df[[time_col]] <= tt)
    idx_t   <- which(df[[time_col]] == tt)
    
    scores_ref   <- run_bod(df[idx_ref, , drop = FALSE])
    pos_t_in_ref <- which(df[idx_ref, time_col] == tt)
    E_S[idx_t] <- scores_ref[pos_t_in_ref]
  }
  
  # 3) intertemporal
  E_I <- run_bod(df)
  
  out <- data.frame(
    id   = df[[id_col]],
    time = df[[time_col]],
    E_C  = E_C,
    E_S  = E_S,
    E_I  = E_I
  )
  
  out <- out[order(out$id, out$time), ]
  
  lag1 <- function(x) c(NA_real_, x[-length(x)])
  
  E_C_lag <- ave(out$E_C, out$id, FUN = lag1)
  E_S_lag <- ave(out$E_S, out$id, FUN = lag1)
  
  # 1. Total Change Index (TC)
  out$TotalChange <- out$E_C / E_C_lag
  
  # 2. Catch-up Index (CU)
  out$CatchUp <- out$E_S / E_S_lag
  
  # 3. Benchmark Shift Index (BS)
  frontier_units <- list()
  for (tt in times) {
    idx <- which(out$time == tt & abs(out$E_C - 1) < 1e-6)
    if (length(idx) > 0) {
      frontier_units[[as.character(tt)]] <- out[idx, ]
    }
  }
  
  out$BenchmarkShift <- NA_real_
  
  for (i in 2:length(times)) {
    t_curr <- times[i]
    t_prev <- times[i - 1]
    
    curr_frontier <- frontier_units[[as.character(t_curr)]]
    prev_frontier <- frontier_units[[as.character(t_prev)]]
    
    if (!is.null(curr_frontier) &&
        !is.null(prev_frontier) &&
        nrow(curr_frontier) > 0 &&
        nrow(prev_frontier) > 0) {
      
      curr_ids <- curr_frontier$id
      prev_ids <- prev_frontier$id
      
      curr_data <- df[df[[time_col]] == t_curr &
                        df[[id_col]] %in% curr_ids,
                      indic_cols, drop = FALSE]
      
      prev_data <- df[df[[time_col]] == t_prev &
                        df[[id_col]] %in% prev_ids,
                      indic_cols, drop = FALSE]
      
      if (nrow(curr_data) > 0 && nrow(prev_data) > 0) {
        curr_means <- colMeans(curr_data, na.rm = TRUE)
        prev_means <- colMeans(prev_data, na.rm = TRUE)
        
        ratios <- curr_means / prev_means
        out$BenchmarkShift[out$time == t_curr] <- mean(ratios, na.rm = TRUE)
      }
    }
  }
  
  out$BenchmarkShift[out$time == times[1]] <- NA_real_
  
  return(out)
}