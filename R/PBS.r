PBS <- function(y,S,B=2000){

	N <- dim(y)[1]
	p <- dim(y)[2]

	reml0 <- REML(y,S)
	e1 <- reml0[[1]][,1]
	tau0 <- reml0[[2]]
	g1 <- tau0^2
	g2 <- 0.5*g1
	
	v1 <- reml0[[1]][,2]^2
	
	G0 <- gmat(g1,g2,p)

	mu0 <- reml0[[1]][,1]

	y.pb <- matrix(numeric(N*p),N)

	T.b <- matrix(numeric(p*B),B)
	
	for(b in 1:B){

		for(i in 1:N){

			yi <- as.vector(y[i,])
			wi <- which(is.na(yi)==FALSE)
			pi <- length(yi)
			mSi <- vmat(S[i,],pi) + G0

			mui <- mu0[wi]
			mSi <- pmat(mSi,wi)

			y.pb[i,] <- ivec2(mvrnorm(1, mui, mSi),wi,p)

		}
		
		remli <- REML(y.pb,S)				# bootstrap estimate
		
		mui <- remli[[1]][,1]
		sei <- remli[[1]][,2]
		taui <- remli[[2]]

		t_new <- mvrnorm(1, mu0, G0)		# pseudo-future observation
		
		T.b[b,] <- (t_new - mui)/sqrt(taui^2+sei^2)
		
		d.b <- ceiling(b/10) - floor(b/10)
		
		if(d.b==0){
		rm1 <- paste0("The ",b,"th bootstrap resampling is completed. (",b,"/",B,")")
		print(rm1)
		}
				
	}

	q1 <- q2 <- numeric(p)
	for(i in 1:p){
		q1[i] <- quantile(T.b[,i],0.025)
		q2[i] <- quantile(T.b[,i],0.975)
	}

	pl <- e1 + q1*sqrt(tau0*tau0+v1)		# prediction interval
	pu <- e1 + q2*sqrt(tau0*tau0+v1)

    R3 <- list("Estimates" = reml0[[1]],"Between-studies_SD" = reml0[[2]],
					  "95%PI" = cbind(pl,pu))
    
	return(R3)
	
}
