pars <- c(iota = x[1],
          rho = rho) #recovery rate

#Parameter penalties

p1x <- ((x[1] > 1.5) * -1.0) * 10^6

outofbounds <- sum(p1x)< 0

if(outofbounds){
  res <- -999999999
} else {
  
  # Catalonia Population: 7675217
  
  susc <- pop - i_0 - r_0
  
  population <- c(susc, i_0, r_0, 0) #Initial condition (population) for ODE
  
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

# .2  0:22, piecewise(15) 0:1, .1, 0:40
interpol <- seq(from = .2, to = .1, length = 17)[-c(1, 17)] 

population <- c(100000 - 300, 300, 0, 0)
pars <- c(iota = 0.2,
          rho = .03)

first_chunk <- 
ode(y = population,
    times = 0:23, func = "derivs", method = "ode45",
    dllname = "SIR", initfunc = "initmod", nout = 0, 
    parms = pars) # Running ODE

for(i in 1:length(interpol)){
  population <- unname(first_chunk[nrow(first_chunk), 2:5])
  pars <- c(iota = interpol[i],
            rho = .03)
  second_chunk <- 
    ode(y = population,
        times = 0:1, func = "derivs", method = "ode45",
        dllname = "SIR", initfunc = "initmod", nout = 0, 
        parms = pars)
  first_chunk <- rbind(first_chunk, second_chunk[2,])
}
population <- unname(first_chunk[nrow(first_chunk), 2:5])
pars <- c(iota = 0.1,
          rho = .03)
third_chunk <- 
  ode(y = population,
      times = 0:39, func = "derivs", method = "ode45",
      dllname = "SIR", initfunc = "initmod", nout = 0, 
      parms = pars)
first_chunk <- rbind(first_chunk, third_chunk[-1,])

ODE_example <- data.frame(date_new = dmy("23/02/2020") + 0:77)

ODE_example <- ODE_example %>% add_column(cases = c(300, diff(first_chunk[, 5])), 
                           popData2020 = 100000, 
                           countriesAndTerritories = "ODE_S3",
                           infected = first_chunk[, 3],
                           recovered = first_chunk[, 4], 
                           dead = 0)
save(ODE_example,  file = "data/second_example_ODE.RData")
