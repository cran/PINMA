# This is a source code file for several computational functions.

library(MASS)

data.edit1 <- function(X1,X2){
  
  r <- is.na(X1[,1])
  X1[r,1] <- 0.001
  X2[r,1] <- 0.01
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    y[,(i-1)] <- log(x1/(x2-x1)) - log(z1/(z2-z1)) 
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    s11 <- 1/z1 + 1/(z2-z1) + 1/x1 + 1/(x2-x1)
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X1[,j]
      w2 <- X2[,j]
      
      s11 <- 1/z1 + 1/(z2-z1) + (x1 - x1) + (w1 - w1)
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X1[,p]
  x2 <- X2[,p]
  
  s11 <- 1/z1 + 1/(z2-z1) + 1/x1 + 1/(x2-x1)
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}

data.edit <- function(study,trt,d,n){

	N <- max(study)
	p <- max(trt)
	
	L <- length(study)

	X1 <- X2 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- trt[i]
  
		X1[k,l] <- d[i]
		X2[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){
  
	if(length(which(X1[i,]==0))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	if(length(which(X1[i,]==X2[i,]))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	}

	X1[is.na(X1[,1]),1] <- 0.001		# data autmentation for the reference arms
	X2[is.na(X2[,1]),1] <- 0.01

	de1 <- data.edit1(X1,X2)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S

	mng2 <- list(y=y,S=S)
  
	return(mng2)
  
}

vmat <- function(q2, p){
  
  i1 <- 1; i2 <- p
  
  Sg <- matrix(numeric(p*p),p,p)
  
  for(i in 1:p){
    
    Sg[i,(i:p)] <- q2[i1:i2]
    
    i1 <- i2 + 1 
    i2 <- i2 + p - i
    
  }
  
  Sg <- Sg + t(Sg); diag(Sg) <- diag(Sg)/2
  
  return(Sg)
  
}

cmat <- function(q2, p){
  
  Sg <- q2*diag(p)
  
  return(Sg)
  
}

pmat <- function(Si, wi){
  
  pl <- length(wi)
  
  R <- matrix(numeric(pl*pl),pl)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[i,j] <- Si[wi[i],wi[j]]
      
    }
  }
  
  return(R)
  
}

imat <- function(Si, wi, p){
  
  pl <- length(wi)
  
  R <- matrix(numeric(p*p),p)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[wi[i],wi[j]] <- Si[i,j]
      
    }
  }
  
  return(R)
  
}

ivec <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- numeric(p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

ivec2 <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- rep(NA,times=p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

gmat <- function(g1,g2,p){
  
  G <- diag(0, p) + g2
  diag(G) <- g1
  return(G)
  
}

QT <- function(x,x0){
  
  x1 <- sort(c(x,x0))
  w1 <- which(x1==as.numeric(x0))
  qt <- 1 - w1/(length(x)+1)
  return(qt)
  
}

fun.I <- function(x){ 
  diag(x)
}
fun.J <- function(x, y = x){ 
  #rep(1, x) %*% t(rep(1, y))
  matrix(1, x, y)
}
fun.e <- function(m, i){
  fun.I(m)[, i]
  #as.numeric((1:m) == i)
}
fun.E <- function(m, i, j){
  fun.e(m, i) %*% t(fun.e(m, j))
}
fun.tilde_P <- function(m){
  (fun.I(m) + fun.J(m)) / 2
}
tr <- function(x){ 
  sum(diag(as.matrix(x)))
  #sum(diag(x))
}
fun.Sum <- function(List){
  rowSums(array(unlist(List), c(dim(List[[1]]), length(List))), dims = 2)
}


