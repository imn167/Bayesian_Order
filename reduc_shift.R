library(nimble)
source('Code_nimble.R')
source('utils.R')
library(tidyverse)
library(coda)


#### Data ####
Delta = c(0, 50,120)
i = 1
sigma = 100
N = 10
sim <- simulated_data(N, Delta[i], sigma, 13)

## fixing the value of the periode 
T1 <- min(sim$Mesure) - 3* sigma
T2 <- max(sim$Mesure) + 3*sigma
periode_T1 = c(periode_T1, T1)
periode_T2 = c(periode_T2, T2)

print(paste0("-------------------------- periode étude [",round(T1),",", round(T2), "]-----------------------------------" ))


raw_data <- read.csv("Reunion_Quinowa/Données_Guillaume.csv", sep = ";", col.names = c("Age_ka", "var"))
raw_data <- raw_data %>% mutate(sd = sqrt(var), depth = 1:17)
N = dim(raw_data)[1]

T1 <- min(raw_data$Age_ka)-3*min(raw_data$sd)
T2 <- max(raw_data$Age_ka)+3*max(raw_data$sd)


  s =runif(1, 0, (T2-T1))
  m1 = runif(1, T1, (T2-s))
  e =  rexp(N-1)
  
  inits = list( s = s , e = e, debut = m1)
  
  
  
  #### Uniform Order ###
  plain_order <- nimbleModel(unif_shift_order, list(N = N), 
                             data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ), 
                             dimensions = list(e = (N-1), mu =N, v = (N-2)),
                             inits = inits)
  cord <-  compileNimble(plain_order)
  conf_ord <-  configureMCMC(plain_order, monitors = c("s", "mu"), thin = 5)
  mcmc <- buildMCMC(conf_ord)
  cmcmc <-  compileNimble(mcmc)

  
  ###changing mcmc ? 
  conf_ord$printSamplers() #Random_Walk samplers for all variables 
  
  samp_order <-  runMCMC(cmcmc, niter = 46000, nburnin = 28000, nchains = 3,
                         progressBar = TRUE, samplesAsCodaMCMC = TRUE)  

  plot(samp_order)  
  
order <- aggregat_samp(samp_order)[, 1:10]

IC_order = data.frame(inf = apply(order, 2, interval_credible)[1, ], 
                      sup = apply(order, 2, interval_credible)[2, ], 
                      estimation = apply(order, 2, median), #MAP meilleur, median voir les proprièté  ? 
                      mesure = sim$Mesure, 
                       names = sim$Names)
IC_order %>% ggplot(aes(x = names, ymin = inf, ymax = sup)) + geom_linerange() +
  geom_point(aes(x = names, y = mesure, colour = "Mesure"))+ theme_imene()





