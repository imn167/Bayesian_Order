library(tidyverse)
library(nimble)
source("Code_nimble.R")
library(ArchaeoPhases)
library(patchwork)
library(ggridges)

#### Creation of a theme (bigger axis.text et some color )

theme_imene <- function() {
  theme_minimal() + 
    theme(axis.text.x = element_text( face = 'bold', color = "#993333", size = 12), 
          axis.text.y = element_text(face = "bold", color = "blue", size = 12), 
          axis.text = element_text(color = "steelblue", face = "bold", size = (10)))
}


#simulation du jeu de données du papier Steier-Rom
simulated_data <- function(n_ages, delta, sigma,inverted = F, n_inv = 2,  half = 1000){
  ### simulation des données 
  set.seed(123) #fixer la souche pour générer des données 
  Y <- rnorm(n_ages, sd = sigma)
  mean <- 1000  + ((1:n_ages) - (N+1) /2) *delta
  M <-  Y + mean
  depth <- seq(200, 400, length.out = n_ages)
  
  ####@
  if (inverted){
    set.seed(123) 
    inverted <- sample(1:n_ages, 2) #inversion de 7 et 8 
    values <- dates[inverted]
    dates[inverted] <- values[n_inv:1]
  }
  ####@
  # dates[5] <- dates[5] + dt/2
  
  dd <- data.frame(Depth = depth, Age = mean, std = rep(sigma, n_ages), Mesure = M,
                   Names = factor(paste0("mu[", 1:n_ages, "]"), levels = paste0("mu[", 1:n_ages, "]")))
  return(dd)
  
}

#### visulisation des densité pour chaque distribution 

########@
comparing <- function(number_age, samp_1, samp_2){
  comparaison <- data.frame(diff = rbind(samp_1[number_age, ], samp_2[number_age, ]), estimator = c('order', 'rademacher')) %>% pivot_longer(!estimator, names_to = 'value', values_to = 'difference')
}

path0 = 'SimuSteierRom/'
Delta = c(0,1,2,5,10,50,100, 120, 150, 180, 200)
N <- 10

#######################################################################################@
plot_data <- list()

dt_var = c()
shift_order = c()

dt_var_rademacher = c()
shift_rademacher = c()

graph_order <- list()
graph_rademacher <- list()

graph_comparing <- list()

periode_T1 = c()
periode_T2 = c()

for (i in seq_along(Delta) ) {
  ###############@
  ### Partie Data
  set.seed(123) #set seed 
  print(Delta[i])
  sim <-simulated_data(N, Delta[i], 100) #donnée avec souche fixer pour Y
  
  
  p <- ggplot(data = sim, mapping = aes(Age, Depth, color = Names)) +geom_point(size = 3, shape = 8) +
    geom_point(aes(Mesure, Depth, color = Names), inherit.aes = F, size = 3) + 
    theme_imene()+ labs(title = paste0("Jeu de données simulé avec delta = ", Delta[i]), subtitle = "Steier_Rom ")
  plot_data[[i]] <- p  
  # view(sim)
  
  T1 <- min(sim$Mesure) - 3* Delta[i]
  T2 <- max(sim$Mesure) + 3*Delta[i]
  periode_T1 = c(periode_T1, T1)
  periode_T2 = c(periode_T2, T2)
  
  print(paste0("-------------------------- periode étude [",T1,",", T2, "]-----------------------------------" ))
  ##############@
  
  ### Ordre classique #####
  plain_order <- nimbleModel(ordre, list(N = N), data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ),
                     inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
  cord <-  compileNimble(plain_order)
  conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
  mcmc <- buildMCMC(conf_ord)
  cmcmc <-  compileNimble(mcmc)

  samp <-  runMCMC(cmcmc, niter = 20000, nburnin = 12000, nchains = 3,
                 progressBar = TRUE, samplesAsCodaMCMC = TRUE)
# plot(samp)
  samp_order <- rbind(samp$chain1, samp$chain2, samp$chain3)
  
  biais_order <- apply(samp_order, 1, function(x) x-sim$Age) / (T2-T1)#distance between bayes estimateur and true age

  apply(samp_order, 2, mean)

  shift_order <- rbind(shift_order, (samp_order[, N] - samp_order[, 1]) / (T2-T1))

  dt <- data.frame(Ic_inf = apply(samp_order, 2, interval_credible)[1,],
                 Ic_sup = apply(samp_order, 2, interval_credible)[2,],
                 estimation = apply(samp_order, 2, mean), true = sim$Age, Mesure = sim$Mesure,
                 Ages = factor(sim$Names, levels = sim$Names) )

  dt_var <- rbind(dt_var, apply(samp_order, 2, var) / (T2-T1)**2 )

  graph <- ggplot(dt) +geom_linerange(aes(x = Ages, ymin = Ic_inf, ymax = Ic_sup)) +
  geom_point(aes(x = Ages, y = estimation), inherit.aes = F) +
  geom_point(aes(x = Ages, y = true, color = 'green'), inherit.aes = F, size = 2) +
    geom_point(aes(x = Ages, y = Mesure, color = 'purple'), inherit.aes = F, size = 2) +
    theme_imene() +
  labs(title = paste0("Ordre Uniforme avec delta = ", Delta[i])) + 
    scale_colour_manual(name = '', 
                        values =c('green'='green','purple'='purple'), labels = c('True Age','Mesured Age'))

  print(graph)
  graph_order[[i]] <- graph
  path = paste0(path0, 'Order/', 'IC_', Delta[i], '.pdf')
  ggsave(path, graph)
  
  
### Bruit Rademacher ####
  
  rademacher_betai <- nimbleModel(beta_pi, list(N=N, alpha = 0, beta = 1), 
                                  data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                                  inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                               b = rbeta(N, 1, N), p = runif(N)),
                                  dimensions = list(v = N))
  
  cord_betai <- compileNimble(rademacher_betai)
  conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "p", "b"), thin = 5)
  mcmc_betai <- buildMCMC(conford_betai)
  cmcmc_betai <- compileNimble(mcmc_betai)
  
  
  samp_betai <- runMCMC(cmcmc_betai, niter = 20000, nburnin = 12000, nchains = 3,
                        progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  pdf(paste0(path0, 'Rademacher_pi/', 'chainplot_', Delta[i], '.pdf'))
  plot(samp_betai)
  dev.off()
  sampBetai <- rbind(samp_betai$chain1[, 21:30], samp_betai$chain2[, 21:30], samp_betai$chain3[, 21:30])
  
  biais_betai <- apply(sampBetai, 1, function(x) x-sim$Age) / (T2-T1) ## distance between bayes estimator and the true age
  
  apply(sampBetai, 2, mean)
  
  shift_rademacher <- rbind(shift_rademacher, (apply(sampBetai, 1, max) - apply(sampBetai, 1, min))/(T2-T1))
  
  
  
  dt <- data.frame(Ic_inf = apply(sampBetai, 2, interval_credible)[1,],
                   Ic_sup = apply(sampBetai, 2, interval_credible)[2,],
                   estimation = apply(sampBetai, 2, mean), true = sim$Age, Mesure = sim$Mesure,
                   Ages = factor(sim$Names, levels = sim$Names) )
  ###variance 
  dt_var_rademacher = rbind(dt_var_rademacher, apply(sampBetai, 2, var) /(T2-T1)**2)
  
  graph <- ggplot(dt) +geom_linerange(aes(x = Ages, ymin = Ic_inf, ymax = Ic_sup)) +
    geom_point(aes(x = Ages, y = estimation), inherit.aes = F) +
    geom_point(aes(x = Ages, y = true, color = 'green'), inherit.aes = F, size = 2) +
    geom_point(aes(x = Ages, y = Mesure, color = 'purple'), inherit.aes = F, size = 2) +theme_imene() +
    labs(title = paste0("Bruit Rademacher avec delta = ", Delta[i])) +
    scale_colour_manual(name = '', 
                        values =c('green'='green','purple'='purple'), labels = c('True Age','Mesured Age'))
  print(graph)
  graph_rademacher[[i]] <- graph
  path = paste0(path0, 'Rademacher_pi/', 'IC_', Delta[i], '.pdf')
  ggsave(path, graph)
  
  
  #######################################@
  ## comparaison des distribtions du biais dans les deux cas 
  comparing_delta <- list()
  for (number_age in 1:N) {
    comparaison <- comparing(number_age, biais_order, biais_betai)
    plot_bias <- ggplot(comparaison, aes(x = difference, group = estimator, color = estimator)) + geom_density() +theme_imene() + 
      labs(title = paste0('Comparaison biais de A', number_age, 'pour delta = ', Delta[i]))
    comparing_delta[[number_age]] <- plot_bias
  }
  graph_comparing[[i]] <- comparing_delta
}

######################################################################################################@@

for (i in c(1,3,5, 7, 9)) {
  print((graph_rademacher[[i]] + graph_rademacher[[i+1]]) / (graph_order[[i]] +graph_order[[i+1]]) )
}

plot_data

dt_var_order<- as.data.frame(dt_var)
dt_var_order = dt_var_order %>% mutate(delta = factor(Delta, levels = Delta))


dt_var_order <- dt_var_order %>% gather(key = 'Age', value = 'Value', -delta) 
dt_var_order <-  dt_var_order %>% mutate(Age = factor(dt_var_order$Age, levels = paste0('mu[',1:N, ']')))

var_order <- ggplot(dt_var_order, aes(x = Age, y = Value, group = delta, color = delta)) +
  geom_line() + theme_imene() + geom_text(aes(label = delta)) 
var_order
ggsave(paste0(path0, 'Order/', 'var.pdf'), var_order)

dt_shift_order <- data.frame(shift_order )
dt_shift_order <- dt_shift_order %>% mutate(delta = factor(Delta, levels = Delta)) %>% 
  pivot_longer(!c(delta), names_to = "sim", values_to = 'value')
ggplot(dt_shift_order, aes(x = value, group = delta, color = delta)) + 
  geom_density() + theme_imene() + xlim(-1,2)

ggplot(dt_shift_order, aes(value, y = delta, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.975), scale = 1, alpha =.3,
    rel_min_height = 1e-8) + theme_minimal() + scale_fill_manual(
      name = "credibility interval at 95%", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0"),
      labels = c("(0, 0.025]","(0.025, 0.975)", "(0.975, 1]")
    ) +theme_imene() + xlim(0,2)

ggsave(paste0(path0, 'Order/', 'shiftdensity.pdf'))

######@

dt_var_Betai<- as.data.frame(dt_var_rademacher)
dt_var_Betai
dt_var_Betai = dt_var_Betai %>% mutate(delta = factor(Delta, levels = Delta)) %>% 
  gather(key = 'Age', value = 'Value', -delta)
dt_var_Betai <-  dt_var_Betai %>% mutate(Age = factor(dt_var_Betai$Age, levels = paste0('mu[',1:N, ']')))
var_betai <- ggplot(dt_var_Betai, aes(x = Age, y = Value, group = delta, color = delta)) +
  geom_line() + theme_imene()  + geom_text(aes(label = delta))
var_betai
ggsave(paste0(path0, 'Rademacher_pi/', 'var.pdf'), var_betai)

dt_shift_rademacher <- data.frame(shift_rademacher)
dt_shift_rademacher <- dt_shift_rademacher %>% mutate(delta = factor(Delta, levels = Delta)) %>% 
  pivot_longer(!delta, names_to = "sim", values_to = 'value')

ggplot(dt_shift_rademacher, aes(x = value, group = delta, color = delta)) +
  geom_density() + theme_imene() + xlim(-1,2)

ggplot(dt_shift_rademacher, aes(value, y = delta, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.975), scale = 1, alpha =.3,
    rel_min_height = 1e-8) + theme_minimal() + scale_fill_manual(
      name = "credibility interval at 95%", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0"),
      labels = c("(0, 0.025]","(0.025, 0.975)", "(0.975, 1]")
    ) + theme_imene()

ggsave(paste0(path0, 'Rademacher_pi/', 'shiftdensity.pdf'))


graph <- graph_comparing[[3]]
graph[[1]] + graph[[N]]
# save(graph_comparing, file = 'difference_graph')
load('difference_graph')
for (i in 1:floor(length(Delta)/2)) {
  print(i)
      print(graph[[i]] + graph[[floor(length(Delta)/2) + i]])      
}

head(dt_shift_rademacher)
dt_shift <- dt_shift_rademacher %>% mutate(value_order = dt_shift_order$value)

ggplot(dt_shift, aes(value, y = delta, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.975), scale = 1, alpha =.3,
    rel_min_height = 1e-8) + 
  stat_density_ridges( aes(value_order, y = delta, fill = 'blue', alpha = .1), inherit.aes = F, 
    geom = "density_ridges_gradient", 
    calc_ecdf = FALSE,
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.975), scale = 1, alpha =.1,
    rel_min_height = 1e-8) +
  theme_imene() + 
  scale_fill_manual(
      name = "Rademacher IC(95%)", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0","#00AAA0" ),
      labels = c("(0, 0.025]","(0.025, 0.975)", "(0.975, 1]", "ordre") 
    ) + theme_imene()

########################################################################################################@@

### Ordre tq Ai = A(i-1) + Z (ui - u(i-1) ) #####

for (delta in c(0,1,2,5,10,50,100, 120, 150, 180, 200)) {
  print(delta)
  sim <-simulated_data(N, delta, 100)
  # view(sim)
  
  T1 <- min(sim$Mesure) - 3* delta -1
  T2 <- max(sim$Mesure) + 3*delta +1 
  
  order_Ai <- nimbleModel(order_A, list(N = N), data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                          inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .7), p = runif(N)), 
                          dimensions = list(v = N))
  
  
  
  cord_Ai <- compileNimble(order_Ai)
  conf_orderAi <- configureMCMC(order_Ai, monitors = c("mu", "Z", "p"), thin = 5)
  mcmc_Ai <- buildMCMC(conf_orderAi)
  cmcmc_Ai <- compileNimble(mcmc_Ai)
  
  
  samp_Ai <- runMCMC(cmcmc_Ai, niter = 20000, nburnin = 12000, nchains = 3,
                     progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  # plot(samp)
  sampAi <- rbind(samp_Ai$chain1[,11:20], samp_Ai$chain2[,11:20], samp_Ai$chain3[,11:20])
  apply(sampAi, 2, mean)
  
  
  dt <- data.frame(Ic_inf = apply(sampAi, 2, interval_credible)[1,],
                   Ic_sup = apply(sampAi, 2, interval_credible)[2,],
                   estimation = apply(sampAi, 2, mean), true = sim$Age,
                   Ages = factor(sim$Names, levels = sim$Names) )
  
  graph <- ggplot(dt) +geom_linerange(aes(x = Ages, ymin = Ic_inf, ymax = Ic_sup)) +
    geom_point(aes(x = Ages, y = estimation), inherit.aes = F) +
    geom_point(aes(x = Ages, y = true), inherit.aes = F, color = 'green', size = 2) +theme_minimal() +
    labs(title = paste0("Ordre Uniforme avec delta = ", delta))
  print(graph)
  path = paste0(path0, 'Order_A/', 'IC_', delta, '.pdf')
  ggsave(path, graph)
  
}





