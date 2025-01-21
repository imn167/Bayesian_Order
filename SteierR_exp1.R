library(nimble)
source("Code_nimble.R")
library(ArchaeoPhases)
library(patchwork)
library(ggridges)
library(tidyverse)
source('utils.R')
library(coda) ### generate trajectory of the mcmc sampling, including the estimated density
#-----------------------------------------------------------------------------@
#Question : - Ajouter la periode d'étude sur le graphique (oui pour les graph des données / mesures)
#           - Utiliser la median, moyenne ?
#           - Tracer le prior pour les lois (posteriori et priori) ? 



#-----------------------------------------------------------------------------@
#### Priors ####
#simulation des loi a prior de Beta Rademacher --> code Python 
path = "./BR_priors/sim_n10.csv"
prior = read.csv(path, row.names = 1)
prior_order = read.csv("./BR_priors/sim_order_n10.csv", row.names = 1)
prior_var <- read.csv("./BR_priors/var_BR_n10.csv", row.names = 1)
prior_var_order <- read.csv("./BR_priors/var_order_n10.csv", row.names = 1)
prior_Et <- read.csv("./BR_priors/Et_prior.csv", row.names = 1)
prior_Et_order <- read.csv("./BR_priors/Et_prior_order.csv", row.names = 1)

#-----------------------------------------------------------------------------@



#####@
path0 = 'SimuSteierRom/'
Delta = c(0,1,2,5,10,50,100,120,150,180,200) #difference between the ages / uniform difference 
N <- 10
sigma <- 100
l <- length(Delta) 

#### Initialisation ####
### probability of having an inversion of 2 consecutives ages :
inversion_prob <- data.frame(prob = pnorm(0, Delta, sqrt(2)*sigma), delta = Delta)
inversion_prob %>% ggplot(aes(delta, prob)) + geom_line() + geom_point() +theme_imene() + ylab("inversion probability")
ggsave(paste0(path0,"inversion_prob.png"))

plot_data <- list() #list of graph for each case of data 

# matrix containing simulation of max(Ai) -min(Ai) for each Delta
shift_order = c() 
shift_rademacher = c()

graph_biais <- list() #IC plot for each delta's value
graph_variance <- list()

graph_comparing <- list() 

graph_post_prior <- list()
graph_post_pi <- list()

periode_T1 = c()
periode_T2 = c()

meanZ <- c()

#loop for each data 
for (i in seq_along(Delta)) {
  print(paste0("------------", Delta[i], '---------------------'))
  #### simulated data #### 
  sim <- simulated_data(N, Delta[i], sigma, 123)
  
  ## fixing the value of the periode 
  T1 <- min(sim$Mesure) - 3* sigma
  T2 <- max(sim$Mesure) + 3*sigma
  periode_T1 = c(periode_T1, T1)
  periode_T2 = c(periode_T2, T2)
  print(paste0("-------------------------- periode étude [",T1,",", T2, "]-----------------------------------" ))
  
  #representation of the likelyhood 
  p <- sim %>% uncount(1000) %>% mutate(value = rnorm(n(), Age, std)) %>% ggplot(aes(x = value, y = Depth, group = Depth)) + 
    geom_density_ridges(scale = 1,fill= pallet[1], alpha = 0.4) + 
    geom_point(data = sim, aes(Mesure, Depth), color = pallet[2], inherit.aes = T, size = 3) + theme_imene() 
  
  plot_data[[i]] <- p
  
  ### initialisation for the order shift 
  s =runif(1, 0, (T2-T1))
  m1 = runif(1, T1, (T2-s))
  e =  rexp(N-1)
  
  inits_order = list( s = s , e = e, debut = m1)
  inits_br = list(s = s, e = e, debut = m1, b = rbeta(N, 1, N), Z1 = rbinom(N, 1, .5), p = runif(N))
  #### Bayesian Model ####
  if (i == 1) {
    #### Uniform Order ###
    plain_order <- nimbleModel(unif_shift_order, list(N = N), 
                               data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ), 
                               dimensions = list(v = (N-2), mu = N, e = (N-1)), inits = inits_order)
    cord <-  compileNimble(plain_order)
    conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
    mcmc <- buildMCMC(conf_ord)
    cmcmc <-  compileNimble(mcmc)
    
    
    #### Beta-Rademacher ###
    rademacher_betai <- nimbleModel(unif_shift_BR_pi, list(N=N, alpha = 0, beta = 1), 
                                    data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                                    inits = inits_br,
                                    dimensions = list(v = (N-2), mu = N, e = (N-1), Z = N, b = N))
    
    cord_betai <- compileNimble(rademacher_betai)
    conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "b", "p"), thin = 5)
    mcmc_betai <- buildMCMC(conford_betai)
    cmcmc_betai <- compileNimble(mcmc_betai)
    
  }
  
  ### setting data 
  cord$setData(list(M = sim$Mesure, tau = sim$std, T1 = T1, T2= T2)) #order
  cord_betai$setData(list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2)) #Beta_pi
  
  #### MCMC run 
  samp_betai <- runMCMC(cmcmc_betai, niter = 30000, nburnin = 22000, nchains = 3,
                        progressBar = TRUE, samplesAsCodaMCMC = TRUE, inits = inits_br)
  samp_order <-  runMCMC(cmcmc, niter = 26000, nburnin = 18000, nchains = 3,
                         progressBar = TRUE, samplesAsCodaMCMC = TRUE, inits = inits_order)
  
  #### MCMC Plot ####
  
  ## looking up the convergence 
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/mu_delta", Delta[i], ".pdf")
  get_mcmc(samp_betai, 2, path_cv)
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/pi_delta", Delta[i], ".pdf")
  get_mcmc(samp_betai, 3, path_cv)
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/Z_delta", Delta[i], ".pdf")
  get_mcmc(samp_betai, 0, path_cv)
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/beta_delta", Delta[i], ".pdf")
  get_mcmc(samp_betai, 1, path_cv)
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/ZB_delta", Delta[i], ".pdf")
  beta <- get_step_br(samp_betai, path_cv, 3)
  path_cv <- paste0(path0, "Comparaison/mcmc_cv/order_delta", Delta[i], ".pdf")
  get_mcmc(samp_order, 0, path_cv)
  
  
  
  ## aggregating chains 
  samp_pi <- aggregat_samp(samp_betai)[, 31:40]
  sampBetai <- aggregat_samp(samp_betai)[, 21:30] #ages from 21-30 
  sampOrder <- aggregat_samp(samp_order)
  sampZ <- aggregat_samp(samp_betai)[, 1:N]
  
  ### mean va latente Z 
  meanZ <- rbind(meanZ, apply(sampZ, 2, mean))
  
  #### biais ####
  biais_order <- apply(sampOrder, 1, function(x) x-sim$Age) / (T2-T1)
  biais_betai <- apply(sampBetai, 1, function(x) x-sim$Age) / (T2-T1)
  
  ic_biais <- data.frame(inf = c(apply(biais_order, 1, interval_credible)[1,], apply(biais_betai, 1, interval_credible)[1,]), 
                         sup = c(apply(biais_order, 1, interval_credible)[2,], apply(biais_betai, 1, interval_credible)[2,]),
                         model = c(rep("order", N), rep("BR", N)), names = rep(sim$Names, 2)
                         )
  graph_biais[[i]] <- ic_biais %>% ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + 
    geom_linerange(position = position_dodge(width = 0.2)) + theme_imene(rotation_x = T)+
    geom_hline(yintercept = 0, linetype = 3)
  
  
  #### shift ####
  shift_order <- rbind(shift_order, (sampOrder[, N] - sampOrder[, 1]) / (T2-T1))
  shift_rademacher <- rbind(shift_rademacher,(sampBetai[, N]-sampBetai[, 1] )/ (T2-T1))
                            #rbind(shift_rademacher, (apply(sampBetai, 1, max) - apply(sampBetai, 1, min))/(T2-T1))
  
  
  
  #### Credibily Interval #### 
  #with estimation and true age 
  IC_order = data.frame(inf = apply(sampOrder, 2, interval_credible)[1, ], 
                        sup = apply(sampOrder, 2, interval_credible)[2, ], 
                        estimation = apply(sampOrder, 2, median), #MAP meilleur, median voir les proprièté  ? 
                        true = sim$Age, 
                        mesure = sim$Mesure, names = sim$Names)
  
  IC_Betai = data.frame(inf = apply(sampBetai, 2, interval_credible)[1, ], 
                        sup = apply(sampBetai, 2, interval_credible)[2, ], 
                        estimation = apply(sampBetai, 2, median), #MAP meilleur, median voir les proprièté  ? 
                        true = sim$Age, 
                        mesure = sim$Mesure, names = sim$Names)
  IC_pi <- data.frame(median = apply(samp_pi, 2, median), names = factor(paste0("pi[", 1:N, "]"), levels = paste0("pi[", 1:N, "]")))
  
  prior_t <- (T2-T1)*prior + T1
  
  ## computing the prior IC(95%) 
  IC_prior <- data.frame(inf = apply(prior_t, 2, interval_credible)[1, ], 
                         sup = apply(prior_t, 2, interval_credible)[2, ], names = sim$Names)
  
  prior_order_t <- (T2-T1) * prior_order + T1
  
  Ic_prior_order <- data.frame(inf = apply(prior_order_t, 2, interval_credible)[1,], 
                               sup = apply(prior_order_t, 2, interval_credible)[2, ], names = sim$Names)
  
  
  #### Comparaison Graph ####
  #using the same scale for the post_prior and _order_BR posterior
  M = max(cbind(Ic_prior_order$sup, IC_prior$sup, IC_Betai$sup))
  m = min(cbind(Ic_prior_order$inf, IC_prior$inf, IC_Betai$inf))
  
  ## comparing the order model and the BR model 
  ic <- (IC_order %>% mutate(model = rep('order', N))) %>% bind_rows(IC_Betai %>% mutate( model = rep("BR", N))) %>% 
    ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + geom_linerange(position = position_dodge(width = 0.2)) + 
    theme_imene(rotation_x =T) + ylim(m, M) + 
    geom_point(aes(x = names, y = mesure, color = "Obs")) + 
    geom_point(aes(x = names, y = true, color = "True")) + labs(color = 'legend') +
    scale_colour_manual(
      values = c("Obs" = pallet[3], "True" = pallet[4], "order" = pallet[5], "BR" = pallet[2])) 
  
  
  
  #### GGplot ####
  ### plotting the BR model posterior and prior 
  ic_post_prior <- IC_Betai %>% mutate(model = rep('BR_posterior', N))  %>% select(inf, sup, names, model) %>% bind_rows(IC_prior %>% mutate(model = "BR_prior")) %>% 
    ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + geom_linerange(position = position_dodge(width = 0.2)) + theme_imene() + ylim(m, M)
  
  ###plotting the pi trajectory (posterior)
  graph_post_pi[[i]] <- IC_pi %>% ggplot(aes(x = names, y = median, group = 1 )) + geom_point() + geom_line() + theme_imene(rotation_x = T) +ylim(.001, .999)
  
  graph_post_prior[[i]] <- ic_post_prior + labs(title = paste0("Delta = ", Delta[i]))
  
  graph_comparing[[i]] <- ic + labs(title = paste0("Delta = ", Delta[i]))
  
  ## variance 
  variance <- data.frame(var_order_post = apply(sampOrder, 2, var) / (T2-T1)^2, 
                         var_BR_post = apply(sampBetai, 2, var)/ (T2-T1)^2, ages = sim$Names)
  
  graph_variance[[i]] <- variance %>% pivot_longer(-ages) %>% ggplot(aes(x = ages, y = value, colour = name, group = name)) + 
    geom_line() +geom_point() +theme_imene(rotation_x = T) +labs(title = paste0("Delta = ", Delta[i]))
  
  
  
  
}#### END OF LOOP #### 

#### Shift Plot ####
ic_shift <- data.frame(inf = c(apply(shift_order, 1, interval_credible)[1,], 
                               apply(shift_rademacher, 1, interval_credible)[1,]), 
                       sup = c(apply(shift_order, 1, interval_credible)[2,], 
                               apply(shift_rademacher, 1, interval_credible)[2,]),
                       model = c(rep("order(uniform_shift)", length(Delta)), rep("BR", length(Delta))), Delta = factor(rep(Delta, 2), 
                      levels = as.character(c(Delta, Delta[length(Delta)] + 1))) 
)
ic_shift <- ic_shift %>% add_row(inf = interval_credible(prior_Et$X0)[1], sup = interval_credible(prior_Et$X0)[2], model = "Et(BR_prior)", Delta = as.factor(201)) %>% 
  add_row(inf = interval_credible(prior_Et_order$X0)[1], sup = interval_credible(prior_Et_order$X0)[2], 
          model = "Et(order_prior)", Delta = as.factor(201))
ic_shift
sapply(ic_shift, class)
graph_shift <-  ic_shift %>% ggplot(aes(x = Delta, ymin = inf, ymax = sup, color = model)) + 
  geom_linerange(position = position_dodge(width = 0.2)) + theme_imene(rotation_x = T)
graph_shift


## the shift here is important in the case where the delta_t is small 
## In fact the shidt that was pint-point for the order 

##----------------------------------------------------------@
meanZ <- as.data.frame(meanZ) %>% mutate(delta = factor(Delta, levels = Delta))
meanZ %>% pivot_longer(!delta) %>%mutate(name = factor(name, levels = paste0("Z[", 1:N, "]"))) %>% 
  ggplot(aes(name, value, group = delta, colour = delta)) + geom_point() + geom_line() + theme_imene()



##------------------------------------------------@
# simulated data visualisation : the distribution of the data + posterior and prior BR IC + Comparaison between BR and order posterior (points mesure and true value)

# have the same scale for mesure (y)
for (i in 1:l) {
  print(plot_data[[i]] +geom_segment(aes(x= periode_T1[i], xend = periode_T2[i], y = sim$Depth[1], yend = sim$Depth[1]), linewidth = 0.5, color = pallet[3]) +
          labs(title = paste0("Delta = ",Delta[i])) + theme(axis.text.x = element_text(angle = 45)) +
    # graph_biais[[i]] + labs(title = paste0("delta = ", Delta[i]))
      graph_post_prior[[i]] + labs(title = paste0("delta = ", Delta[i])) + theme(axis.text.x = element_text(angle = 45)) +
      graph_comparing[[i]] + labs(title = paste0("delta = ", Delta[i])) + theme(axis.text.x = element_text(angle = 45))  )
  #Biais graph doesn't add any new information 
  ggsave(paste0(path0, "Comparaison/Delta_", Delta[i], ".png"), width= 13, height = 7)
}


for (i in 1:l) {
  print(graph_comparing[[i]])
}
### Variance a prior pour BR +++ importante 
##--------------------------------------------------------@
var_prior <- data.frame( var_order_prior = prior_var_order$X0, var_BR_prior = prior_var$X0, ages = variance$ages) %>% pivot_longer(!ages)
var_prior %>% ggplot(aes(ages, value, group = name, colour= name)) + geom_point() + geom_line() + theme_imene() + theme(axis.text.x = element_text(angle = 45))

for (i in 1:l) {
  print(graph_variance[[i]] + graph_post_pi[[i]] )
  ggsave(paste0(path0, "Comparaison/Var_pi_Delta_", Delta[i], ".png"), width= 13, height = 7)
}



#---------------------------------------------------------------@
###### Testing 

nimble_sort(runif(12))
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------@


