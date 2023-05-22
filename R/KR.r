KR <- function(y, S){
  
  result_reml <- REML(y,S)
  tau_reml <- result_reml[[2]]
  kr <- KR0(y = y, S = S, tau2 = tau_reml^2)
  
  kr[["Expected"]][["df"]][kr[["Expected"]][["df"]] < 3] <- 3		# truncation of the DF
  
  KR_E <- list("Coefficients" = cbind("Coef." = result_reml[[1]][, 1], 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "95%CU" = result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml)
  
  #
  
  pl <- KR_E[[1]][,1] - qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2)
  pu <- KR_E[[1]][,1] + qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2)

  R3 <- list("Estimates" = cbind("Coef." = result_reml[[1]][, 1], 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "95%CU" = result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml,
					  "95%PI" = cbind(pl,pu))
    
  return(R3)
  
}

