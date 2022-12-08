ll_ode <- function(x, # parameter vector
                   y, # data
                   devs, #standard deviation
                   pop, #population of the country
                   rho,
                   i_0,# Maximum infected people
                   r_0){ # Maximum recovered people
  
  pars <- c(iota = x[1],
            rho = rho) #recovery rate

  #Parameter penalties
  
  p1x <- ((x[1] > 1.5) * -1.0) * 10^6
  p2x <- ((x[2] > 1.5) * -1.0) * 10^6
  ptimex <- ((round(x[3]) > (nrow(y) - 1)) * -1.0) * 10^6
  ptimen <- ((round(x[3]) < 2) * -1.0) * 10^6
  
  outofbounds <- sum(p1x, p2x, ptimex, ptimen) < 0
  
  if(outofbounds){
    res <- -999999999
  } else {
    
    # Catalonia Population: 7675217
    
    susc <- pop - i_0 - r_0
    
    population <- c(susc, i_0, r_0, 0) #Initial condition (population) for ODE
    
    z1 <- ode(y = population,
             times = 0:round(x[3]), func = "derivs", method = "ode45",
             dllname = "SIR", initfunc = "initmod", nout = 0, 
             parms = pars) # Running ODE
    
    # Changing values for the break
    
    pars <- c(iota = x[2],
              rho = rho)
    population <- z1[(round(x[3]) + 1), -1]
    z2 <- ode(y = population,
              times = round(x[3]):nrow(y), func = "derivs", method = "ode45",
              dllname = "SIR", initfunc = "initmod", nout = 0, 
              parms = pars) # Running ODE
    
    colnames(z1)[2:5] <- c("S", "I", "R", "Pos_ac")
    colnames(z2)[2:5] <- c("S", "I", "R", "Pos_ac")
    
    z <- rbind(z1, z2[-1, ])
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    Pos_ac <- (y$cases)
    
    res <- #This is the log-likelihood
      sum(dnorm(diff(c(0, z$Pos_ac)) - Pos_ac, sd = devs, log = T))
    
    rm(z)
    rm(y)
    
  }
  res
}
