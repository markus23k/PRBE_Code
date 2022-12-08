fitting_procedure <- function(seeds, fun, data, stdev, pop, rho, i_0, r_0){
  
  best <- -999999999 #pre-alocating LL
  
  sols <- NA #Pre-alocating the good (at a distance of 2 LL from the best) solutions
  
  sims <- ncol(seeds)
  condition <- T 
  # This condition stops the fitting procedure if the solutions don't get 2 LL
  # points better and there are some good solutions
  
  round <- 1
  
  while(condition){
    parall <- list()
    for(k in 1:ncol(seeds)){
      
      fit <- optim(par = seeds[, k], fn = fun, y = data,  
                   devs = stdev, pop = pop, rho = rho, i_0 = i_0, r_0 = r_0,
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
    
    
    pc3 <- sols > 10 # Maybe Change condition pc3
    
    
    condition <- pc1 * pc2 * pc3
    condition <- !condition
    
    # We select the good (or not so bad solutions for the next round)
    if(sols < 30){
      index <- order(logl, decreasing = T)[1:25]
    } else {
      index <- order(logl, decreasing = T)[1:sols]
    }
    
    last_round_params <- matrix(NA, nrow = sims, ncol = nrow(seeds))
    n <- 1
    parmat <- matrix(NA, nrow = length(index), ncol = nrow(seeds))
    for(i in 1:sims){
      last_round_params[n, ] <- lhs[[i]]$par
      n <- n + 1
    }
    
    parmat <- last_round_params[index, ]
    
    # I sample the parameter combination independently. Thus, I am breaking
    # correlations between parameters.
    
    seeds <- t(apply(X = last_round_params, MARGIN = 2, sample, size = sims, replace = T))
    # print(seeds)
    # print(parmat)
    # print(last_round_params)
    # # last_round_params <- parmat
    
    rm(parmat)
    rm(lhs)
    
    # print(paste0("Round: ", round, ";  Best likelihood: ", best2, 
    #              "; Solutions: ", sols, "\n"))
    # 
    best <- best2
    round <- round + 1
    
    if (round > 10) break
  } 
  last_round_params[index, ]
}
