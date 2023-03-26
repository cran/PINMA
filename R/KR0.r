KR0 <- function(y, S, tau2){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  matrix.indicator <- !(is.na(y))
  
  listlist.y_X_Siinv_tildeP <- lapply(1:N, function(i){
    indicator_i <- matrix.indicator[i, ]
    y_i <- y[i, indicator_i]
    X_i <- fun.I(p)[indicator_i, , drop = FALSE]
    c_i <- length(y_i)
    S_i <- matrix(NA, c_i, c_i)
    S_i[lower.tri(S_i, diag = TRUE)] <- c(na.omit(S[i, ]))
    lower_S_i <- S_i
    S_i <- t(S_i)
    S_i[lower.tri(S_i)] <- lower_S_i[lower.tri(lower_S_i)]
    tilde_P_i <- fun.tilde_P(c_i)
    #Siinv_i <- solve(tau2 * tilde_P_i + S_i)
    Siinv_i <- ginv(tau2 * tilde_P_i + S_i)
    return(list(y_i, X_i, Siinv_i, tilde_P_i))
  })
  
   Sum_xppx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_yppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[1]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xpx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xpy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[2]]
  }))
  Sum_xy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[1]]
  }))
  
  # Phi, P, Q
  #Phi <- solve(Sum_xx)
  Phi <- ginv(Sum_xx)
  P <- - Sum_xpx
  Q <- Sum_xppx
  
  #I_E
  I_E_1 <- sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      tr(z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]])
    })
  ))
  I_E_2 <- tr(2 * (Phi %*% Q) - Phi %*% P %*% Phi %*% P)
  I_E <- (1 / 2) * (I_E_1 - I_E_2)
  
  # I_O
  I_O_1 <- - sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      sum((z[[3]] %*% z[[4]] %*% z[[3]]) * z[[4]])
    })
  ))
  I_O_2 <- sum(((- Phi %*% P %*% Phi) * P) + (Phi * (2 * Q)))
  I_O_3 <- 2 * c(
    Sum_yppy - t(Sum_xpy) %*% Phi %*% Sum_xpy - 2 * t(Sum_xy) %*% Phi %*% Sum_xppy + 2 * t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpy + t(Sum_xy) %*% Phi %*% Sum_xppx %*% Phi %*% Sum_xy - t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xy
  )
  I_O <- (1 / 2) * (I_O_1 + I_O_2 + I_O_3)
  
  # Phi_A_E, Phi_A_O
  Phi_A_E <- Phi + 2 * Phi %*% ((1 / I_E) * (Q - P %*% Phi %*% P)) %*% Phi
  Phi_A_O <- Phi + 2 * Phi %*% ((1 / I_O) * (Q - P %*% Phi %*% P)) %*% Phi
  
  vec.W <- 1 / c(I_E, I_O)
  fun.m_lambda <- function(j){
    L_j <- as.matrix(fun.e(p, j))
    l_j <- 1
    #Theta_j <- L_j %*% solve(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    Theta_j <- L_j %*% ginv(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    vec.A1_j <- vec.W * (tr(Theta_j %*% Phi %*% P %*% Phi))^2
    vec.A2_j <- vec.W * tr(Theta_j %*% Phi %*% P %*% Phi %*% Theta_j %*% Phi %*% P %*% Phi)
    vec.g_j <- ((l_j + 1) * vec.A1_j - (l_j + 4) * vec.A2_j) / ((l_j + 2) * vec.A2_j)
    vec.c1_j <- vec.g_j / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c2_j <- (l_j - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c3_j <- (l_j + 2 - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.Estar_j <- 1 / (1 - vec.A2_j / l_j)
    vec.B_j <- (vec.A1_j + 6 * vec.A2_j) / (2 * l_j)
    vec.Vstar_j <- (2 / l_j) * (1 + vec.c1_j * vec.B_j) / ((1 - vec.c2_j * vec.B_j)^2 * (1 - vec.c3_j * vec.B_j))
    vec.rho_j <- vec.Vstar_j / (2 * vec.Estar_j^2)
    vec.m_j <- 4 + (l_j + 2) / (l_j * vec.rho_j - 1)
    vec.lambda_j <- vec.m_j / (vec.Estar_j * (vec.m_j - 2))
    return(c(m_E = vec.m_j[1], lambda_E = vec.lambda_j[1], m_O = vec.m_j[2], lambda_O = vec.lambda_j[2]))
  }
  matrix.m_lambda <- sapply(1:p, fun.m_lambda)
  
  SE_E <- sqrt(pmax(0, diag(Phi_A_E)))
  SE_O <- sqrt(pmax(0, diag(Phi_A_O)))
  df_E <- matrix.m_lambda[1, ]
  lambda_E <- matrix.m_lambda[2, ]
  df_O <- matrix.m_lambda[3, ]
  lambda_O <- matrix.m_lambda[4, ]
  
   if(any(diag(Phi) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi` is negative. ")
  }
  if(any(diag(Phi_A_E) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_E` is negative. ")
  }
  if(any(diag(Phi_A_O) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_O` is negative. ")
  }
  if(any(df_E < 0)){
    warning("At least one element of the vector `df_E` is negative. ")
  }
  if(any(df_O < 0)){
    warning("At least one element of the vector `df_O` is negative. ")
  }
  
  #
  return(
    list(Expected = list(SE = SE_E, df = df_E, lambda = lambda_E), 
         Observed = list(SE = SE_O, df = df_O, lambda = lambda_O))
  )
  
}


