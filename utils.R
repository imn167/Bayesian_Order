library(tidyverse)
library(ggridges)
library(patchwork)


#### Nimble options 
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

#### creation theme pour ggplot

theme_imene <- function(rotation_x = F) {
  tt <-  theme_minimal() + 
    theme(axis.text.x = element_text( face = 'bold', color = "#993355", size = 12), 
          axis.text.y = element_text(face = "bold", color = "#993355", size = 12),
          axis.title.x = element_text(face = 'bold', size = 14, color = "black"),
          axis.title.y = element_text(face = 'bold', size = 14, color = "black"))
  if (rotation_x) {
    tt <- theme_minimal() + 
      theme(axis.text.x = element_text( face = 'bold', color = "#993355", size = 12, angle = 90), 
            axis.text.y = element_text(face = "bold", color = "#993355", size = 12),
            axis.title.x = element_text(face = 'bold', size = 14, color = "black"),
            axis.title.y = element_text(face = 'bold', size = 14, color = "black"))
  }
  return(tt)
}


simulated_data <- function(n_ages, delta, sigma, seed,  half = 1000){
  ### simulation des données 
  set.seed(seed) #fixer la souche pour générer des données 
  Y <- rnorm(n_ages, sd = sigma)
  mean <- half  + ((1:n_ages) - (n_ages+1) /2) *delta
  M <-  Y + mean
  depth <- seq(200, 400, length.out = n_ages)
  
  #######################@
  dd <- data.frame(Depth = depth, Age = mean, std = rep(sigma, n_ages), Mesure = M,
                   Names = factor(paste0("mu[", 1:n_ages, "]"), levels = paste0("mu[", 1:n_ages, "]")))
  return(dd)
  
}

##############------------------------------------------------------------------------------------------@
#### COmputing bayesian simulation with Nimble 

#init function 
init <- function(N) {
  list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
       b = rbeta(N, 1, N), p = runif(N))
}





##############------------------------------------------------------------------------------------------@
### color pallet 
pallet <- c("#FFF0AAAA", "#0000FFA0", "#00AAA0", "#D44D20", "#9DDF3E", "#3BBFDF", "#F3EC5E", "#ED5524")


##############------------------------------------------------------------------------------------------@


### aggregate the MARKOV CHAIN 

aggregat_samp <- function(samp) {
  n_chain = length(samp)
  samp_chain = c()
  for (i in 1:n_chain) {
    samp_chain = rbind(samp_chain, samp[[i]])
    
  }
  return(samp_chain)
} ##return the samp_chain after MCMC algo 

### element to return from simulation 

samp_properties <- function(samp_chain, data) { 
  ### construction the bias 
  biais = apply(samp_chain, 1, function(x) x-data[, 3]) 
  
}


##-------------------------------------------------------@
#mcmc convergence 
get_mcmc <- function(samp, k, path) {
  m = length(samp) #m chains with diff init
  L <- list() #list de mcmc objects
  c1 <- k*10 +1
  c2 <- (k+1)*10
  
  
  ## Gelman-Rubin : Rhat
  Rhat <- c()
  
  for (i in 1:m) {
    L[[i]] <- coda::mcmc(samp[[i]][, c1: c2])
    
  }
  mcmc_list <- mcmc.list(L)
  pdf(path)
  plot(mcmc_list)
  autocorr.plot(mcmc_list)
  gelman.plot(mcmc_list, autoburnin = F)
  dev.off()
  
  
  
}


#-------------------------------------------------------@
#### Visualisation ####






#### config matrix of data 

config_matrix <- function(mat, factor_vector, var_names) {
  dt <- data.frame(mat)
  dt <- dt %>% mutate(group = factor(factor_vector, levels = factor_vector))
  dt <- dt %>% pivot_longer(!group, names_to = 'Age', values_to = 'value') %>% 
    mutate(Age = factor(dt$Age, levels = paste0('mu[',1:N, ']')))
  
}


























































