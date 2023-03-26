REML <- function(y,S,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  mu <- rnorm(p)	# initial values
  g1 <- 0.2
  g2 <- 0.1
  
  Qc0 <- c(mu,g1,g2)
  
  LL1 <- function(g){
    
    #G <- gmat(g,g2,p)
    G <- gmat(g,(g/2),p)
    
    ll1 <- 0; XWX <- gmat(0,0,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      mui <- mu[wi]
      Gi <- pmat(G,wi)
      
      B1 <- (yi - mui)
      B2 <- ginv(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      A3 <- pl * log(2*pi)
      
      XWX <- XWX + imat(B2,wi, p)
      
      ll1 <- ll1 + A1 + A2 + A3
      
    }
    
    ll2 <- ll1 + log(det(XWX)) - p*log(2*pi)
    
    return(ll2)
    
  }
  
  LL2 <- function(g){
    
    G <- gmat(g1,g,p)
    
    ll1 <- 0; XWX <- gmat(0,0,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      mui <- mu[wi]
      Gi <- pmat(G,wi)
      
      B1 <- (yi - mui)
      B2 <- ginv(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      A3 <- pl * log(2*pi)
      
      XWX <- XWX + imat(B2,wi, p)
      
      ll1 <- ll1 + A1 + A2 + A3
      
    }
    
    ll2 <- ll1 + log(det(XWX)) - p*log(2*pi)
    
    return(ll2)
    
  }
  
  for(itr in 1:maxitr){
    
    A1 <- numeric(p)
    A2 <- matrix(numeric(p*p),p)
    
    G <- gmat(g1,g2,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- ginv(Gi + Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
    mu <- A1 %*% ginv(A2)
    g1 <- optimize(LL1, lower = 0, upper = 5)$minimum
    g2 <- 0.5*g1
    
    V.mu <- ginv(A2)
    
    Qc <- c(mu,g1,g2)
    
    rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
    if(max(rb) < 10^-4) break
    
    Qc0 <- Qc
    
  }
  
  
  SE <- sqrt(diag(V.mu))
  
  R1 <- as.vector(mu)
  R2 <- as.vector(SE)
  R3 <- as.vector(mu - qnorm(.975)*SE)
  R4 <- as.vector(mu + qnorm(.975)*SE)
  
  R5 <- cbind(R1,R2,R3,R4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU")
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
  
  #return(R8)
  return(
    list("Coefficients"=R5,"Between-studies_SD"=R6)
  )
  
}


