ll_ode <- function(x, # parameter vector
                   y, # data
                   devs, #standard deviation
                   pop, #population of the country
                   i_lim,# Maximum infected people
                   r_lim){ # Maximum recovered people
  
  pars <- c(iota = x[1], #per capita infectivity
            rho = x[2]) #recovery rate
  
  #Parameter penalties
  penaltiesn <- ((x > 1e-5) * 1.0  - 1) * 10^6
  
  p1x <- ((x[1] > 1.5) * -1.0) * 10^6
  p2x <- ((x[2] > 0.5) * -1.0) * 10^6
  p2n <- ((x[2] < 0.005) * -1.0) * 10^6
  p3x <- ((x[3] > i_lim) * -1.0) * 10^6 
  p4x <- ((x[4] > r_lim) * -1.0) * 10^6 
  
  outofbounds <- sum(penaltiesn + p1x + p2x + p2n + p3x + p4x)< 0
  
  if(outofbounds){
    res <- -999999999
  } else {
    
    # Catalonia Population: 7675217
    
    susc <- pop - x[3] - x[4]
    
    population <- c(susc, x[3], x[4], 0) #Initial condition (population) for ODE
    
    z <- ode(y = population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             dllname = "SIR", initfunc = "initmod", nout = 0, 
             parms = pars) # Running ODE
    
    colnames(z)[2:5] <- c("S", "I", "R", "Pos_ac")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    Pos_ac <- (y$cases)
    
    res <- #This is the log-likelihood
      sum(dnorm(diff(c(0, z$Pos_ac)) - Pos_ac, sd = devs, log = T))
    
    rm(z)
    rm(y)
    
    res}
}