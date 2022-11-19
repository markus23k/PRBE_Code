fitting_procedure <- function(seeds, fun, data, stdev, pop, i_lim, r_lim){
  
  best <- -999999999 #pre-alocating LL
  
  sols <- NA #Pre-alocating the good (at a distance of 2 LL from the best) solutions
  
  condition <- T 
  # This condition stops the fitting procedure if the solutions don't get 2 LL
  # points better and there are some good solutions
  
  round <- 1
  
  while(condition){
    parall <- list()
    for(k in 1:ncol(seeds)){
      
      fit <- optim(par = seeds[, k], fn = fun, y = data,  
                   devs = stdev, pop = pop, i_lim = i_lim, r_lim = r_lim,
                   control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
      
      
      parall[[k]] <- fit
      rm(fit)
    }
    
    lhs <- parall
    
    rm(parall) #Avoiding memory loss
    
    # Recovering LL 
    logl <- rep(NA, sims)
    for(i in 1:sims) logl[i] <- lhs[[i]]$value
    
    # Evaluating conditions to stop the loop
    best2 <- max(logl, na.rm = T)
    pc1 <- best < best2
    pc2 <- best > (best2 - 2)
    
    sols <- sum(logl > (max(logl, na.rm = T) - 2), na.rm = T)
    pc3 <- sols > 10
    
    
    condition <- pc1 * pc2 * pc3
    condition <- !condition
    
    # We select the good (or not so bad solutions for the next round)
    if(sols < 10){
      index <- order(logl, decreasing = T)[1:25]
    } else {
      index <- order(logl, decreasing = T)[1:sols]
    }
    
    n <- 1
    parmat <- matrix(NA, nrow = length(index), ncol = nrow(seeds))
    for(i in index){
      parmat[n, ] <- lhs[[i]]$par
      n <- n + 1
    }
    
    # I sample the parameter combination independently. Thus, I am breaking
    # correlations between parameters.
    
    seeds <- t(apply(X = parmat, MARGIN = 2, sample, size = sims, replace = T))
    
    last_round_params <- parmat
    
    rm(parmat)
    rm(lhs)
    
    # print(paste0("Round: ", round, ";  Best likelihood: ", best2, 
    #              "; Solutions: ", sols, "\n"))
    # 
    best <- best2
    round <- round + 1
    
    if (round > 10) break
  } 
  last_round_params
}
