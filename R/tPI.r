tPI <- function(y, S){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  A8 <- numeric(p)

  for(i in 1:N){

	vi <- diag(vmat(S[i,],p))
	A8[which(vi<1000)] <- A8[which(vi<1000)] + 1

  }
  
  reml0 <- REML(y,S)

  e1 <- reml0[[1]][,1]
  v1 <- reml0[[1]][,2]^2
  tau2 <- reml0[[2]]^2

  pl <- e1 - qt(0.975,df=N-A8-1)*sqrt(tau2+v1)		# prediction interval
  pu <- e1 + qt(0.975,df=N-A8-1)*sqrt(tau2+v1)

  R3 <- list("Estimates" = reml0[[1]],"Between-studies_SD" = reml0[[2]],
					  "95%PI" = cbind(pl,pu))
    
  return(R3)
  
}



