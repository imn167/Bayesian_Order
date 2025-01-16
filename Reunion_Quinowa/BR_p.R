#################################################################################@
###### In this script we test the model BR with a p fixed p = .5
###### - Test on simulated data (similar to high resolution delta_t small )
##### - Test on the dataset Guillaume sended 

setwd("/Users/E204701B/Desktop/Testing_Priors/Bayesian_Order/Reunion_Quinowa") 

# NOTE: All packages needed are either in the Code_nimble.R file or the utils.R 
#################################################################################@
source("../Code_nimble.R")
source("../utils.R")



##==========================================================================@
#### Priors ####
path = "../BR_priors/BR0.5_sim_n10.csv"
prior_beta = read.csv(path, row.names = 1)
prior_order = read.csv("../BR_priors/sim_order_n10.csv", row.names = 1)
# prior_var <- read.csv("../BR_priors/var_BR_n10.csv", row.names = 1)
# prior_var_order <- read.csv("../BR_priors/var_order_n10.csv", row.names = 1)
# prior_Et <- read.csv("../BR_priors/Et_prior.csv", row.names = 1)
# prior_Et_order <- read.csv("../BR_priors/Et_prior_order.csv", row.names = 1)


##### List #####
N= 10
sigma = 100
Delta = c(0,1,2,5,10,50,100, 120, 150, 180, 200)

plot_data <- list() #list of graph for each case of data 

# matrix containing simulation of max(Ai) -min(Ai) for each Delta
shift_order = c() 
shift_rademacher = c()

graph_comparing <- list() 
graph_variance <- list()
graph_post_prior <- list()

meanZ = c()

periode_T1 = c()
periode_T2 = c()


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
    geom_point(data = sim, aes(Mesure, Depth), color = pallet[2], inherit.aes = T, size = 3) + theme_imene() +
    geom_segment(aes(x = T1, xend = T2, y = sim$Depth[1], yend = sim$Depth[1]), inherit.aes = F, color = pallet[3])
    
  
  plot_data[[i]] <- p
  
  #The number of samples doesn't change, we wil only call the model once and then setdata for the followinf cases 
  #### Bayesian Model ####
  if (i == 1) {
    
    #### Uniform Order ###
    plain_order <- nimbleModel(ordre, list(N = N), data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ),
                               inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
    cord <-  compileNimble(plain_order)
    conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
    mcmc <- buildMCMC(conf_ord)
    cmcmc <-  compileNimble(mcmc)
    
    
    #### Beta-Rademacher ###
    rademacher_beta <- nimbleModel(inv_beta, list(N=N, p= .5), 
                                    data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                                    inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                                 b = rbeta(N, 1, N)),
                                    dimensions = list(v = N))
    
    cord_beta <- compileNimble(rademacher_beta)
    conford_beta <- configureMCMC(rademacher_beta, monitors = c("mu", "Z", "b"), thin = 5)
    mcmc_beta <- buildMCMC(conford_beta)
    cmcmc_beta <- compileNimble(mcmc_beta)
    
    samp_order <-  runMCMC(cmcmc, niter = 26000, nburnin = 18000, nchains = 3,
                           progressBar = TRUE, samplesAsCodaMCMC = TRUE)
    plot(samp_order)
    
  }
  
  ### setting data 
  cord$setData(list(M = sim$Mesure, tau = sim$std, T1 = T1, T2= T2)) #order
  cord_beta$setData(list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2)) #Beta_pi
  
  
  #### MCMC run 
  samp_beta <- runMCMC(cmcmc_beta, niter = 30000, nburnin = 22000, nchains = 3,
                        progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  samp_order <-  runMCMC(cmcmc, niter = 26000, nburnin = 18000, nchains = 3,
                         progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  #### MCMC Plot ####
  
  ## looking up the convergence 
  path_cv <- paste0("mcmc_plot/Z_delta", Delta[i], ".pdf")
  get_mcmc(samp_beta, 0, path_cv)
  path_cv <- paste0("mcmc_plot/beta_delta", Delta[i], ".pdf")
  get_mcmc(samp_beta, 1, path_cv)
  path_cv <- paste0("mcmc_plot/mu_delta", Delta[i], ".pdf")
  get_mcmc(samp_beta, 2, path_cv)
  path_cv <- paste0("mcmc_plot/ZB_delta", Delta[i], ".pdf")
  beta <- get_step_br(samp_beta, path_cv, 3)
  
  
  ## aggregating chains 
  sampBeta <- aggregat_samp(samp_beta)[, (2*N+1):(3*N)] 
  sampZ <- aggregat_samp(samp_beta)[, 1:N]
  sampOrder <- aggregat_samp(samp_order)
  

  ### mean Z 
  meanZ = rbind(meanZ, apply(sampZ, 2, mean))
  
  #### shift ####
  shift_order <- rbind(shift_order, (sampOrder[, N] - sampOrder[, 1]) / (T2-T1))
  # shift_rademacher <- rbind(shift_rademacher, (apply(sampBeta, 1, max) - apply(sampBeta, 1, min))/(T2-T1))
  shift_rademacher <- rbind(shift_rademacher, (sampBeta[, N]- sampBeta[, 1]) / (T2-T1))
  
  
  
  #### Credibily Interval #### 
  #with estimation and true age 
  IC_order = data.frame(inf = apply(sampOrder, 2, interval_credible)[1, ], 
                        sup = apply(sampOrder, 2, interval_credible)[2, ], 
                        estimation = apply(sampOrder, 2, median), #MAP meilleur, median voir les proprièté  ? 
                        true = sim$Age, 
                        mesure = sim$Mesure, names = sim$Names)
  
  IC_Beta = data.frame(inf = apply(sampBeta, 2, interval_credible)[1, ], 
                        sup = apply(sampBeta, 2, interval_credible)[2, ], 
                        estimation = apply(sampBeta, 2, median), #MAP meilleur, median voir les proprièté  ? 
                        true = sim$Age, 
                        mesure = sim$Mesure, names = sim$Names)
  
  
  prior_t <- (T2-T1)*prior_beta + T1
  
  ## computing the prior IC(95%) 
  IC_prior <- data.frame(inf = apply(prior_t, 2, interval_credible)[1, ], 
                         sup = apply(prior_t, 2, interval_credible)[2, ], names = sim$Names)
  
  prior_order_t <- (T2-T1) * prior_order + T1
  
  Ic_prior_order <- data.frame(inf = apply(prior_order_t, 2, interval_credible)[1,], 
                               sup = apply(prior_order_t, 2, interval_credible)[2, ], names = sim$Names)
  
  
  #### Comparaison Graph ####
  #using the same scale for the post_prior and _order_BR posterior
  M = max(cbind(Ic_prior_order$sup, IC_prior$sup, IC_Beta$sup))
  m = min(cbind(Ic_prior_order$inf, IC_prior$inf, IC_Beta$inf))
  
  ## comparing the order model and the BR model 
  ic <- (IC_order %>% mutate(model = rep('order', N))) %>% bind_rows(IC_Beta %>% mutate( model = rep("BR", N))) %>% 
    ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + geom_linerange(position = position_dodge(width = 0.2)) + 
    theme_imene(rotation_x =T) + ylim(m, M) + 
    geom_point(aes(x = names, y = mesure, color = "Obs")) + 
    geom_point(aes(x = names, y = true, color = "True")) + labs(color = 'legend') +
    scale_colour_manual(
      values = c("Obs" = pallet[3], "True" = pallet[4], "order" = pallet[5], "BR" = pallet[2])) 
  
  #### variance ####
  variance <- data.frame(var_order_post = apply(sampOrder, 2, var) / (T2-T1)^2, 
                         var_BR_post = apply(sampBeta, 2, var)/ (T2-T1)^2, ages = sim$Names)
  
  graph_variance[[i]] <- variance %>% pivot_longer(-ages) %>% ggplot(aes(x = ages, y = value, colour = name, group = name)) + 
    geom_line() +geom_point() +theme_imene(rotation_x = T) +labs(title = paste0("Delta = ", Delta[i]))
  
  
  #### GGplot ####
  ### plotting the BR model posterior and prior 
  ic_post_prior <- IC_Beta %>% mutate(model = rep('BR_posterior', N))  %>% select(inf, sup, names, model) %>% 
    bind_rows(IC_prior %>% mutate(model = "BR_prior")) %>% 
    ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + 
    geom_linerange(position = position_dodge(width = 0.2)) + theme_imene(rotation_x = T) + ylim(m, M)
  
  ###plotting the pi trajectory (posterior)
  graph_post_prior[[i]] <- ic_post_prior + labs(title = paste0("Delta = ", Delta[i]))
  
  graph_comparing[[i]] <- ic + labs(title = paste0("Delta = ", Delta[i]))
  
  
}#### END OF LOOP ####


ic_shift <- data.frame(inf = c(apply(shift_order, 1, interval_credible)[1,], 
                               apply(shift_rademacher, 1, interval_credible)[1,]), 
                       sup = c(apply(shift_order, 1, interval_credible)[2,], 
                               apply(shift_rademacher, 1, interval_credible)[2,]),
                       model = c(rep("order", length(Delta)), rep("BR", length(Delta))), Delta = factor(rep(Delta, 2), 
                                                                                                        levels = as.character(c(Delta, Delta[length(Delta)] + 1))) 
)
ic_shift <- ic_shift %>% add_row(inf = interval_credible(prior_Et$X0)[1], sup = interval_credible(prior_Et$X0)[2], model = "Et(BR_prior)", Delta = as.factor(201)) %>% 
  add_row(inf = interval_credible(prior_Et_order$X0)[1], sup = interval_credible(prior_Et_order$X0)[2], model = "Et(order_prior)", Delta = as.factor(201))
ic_shift
sapply(ic_shift, class)
graph_shift <-  ic_shift %>% ggplot(aes(x = Delta, ymin = inf, ymax = sup, color = model)) + 
  geom_linerange(position = position_dodge(width = 0.2)) + theme_imene(rotation_x = T)
graph_shift


for (i in seq_along(Delta)) {
print(graph_comparing[[i]] + graph_post_prior[[i]])
ggsave(paste0("expA_plot/fixed_delta_", Delta[i], ".png"), width= 13, height = 7)
  
}

meanZ <- as.data.frame(meanZ) %>% mutate(delta = factor(Delta, levels = Delta))
meanZ %>% pivot_longer(!delta) %>%mutate(name = factor(name, levels = paste0("Z[", 1:N, "]"))) %>% 
  ggplot(aes(name, value, group = delta, colour = delta)) + geom_point() + geom_line() + theme_imene()


plot_data[[11]] + data.frame(mean = apply(samp_zi, 2, mean), names = paste0("Z[", 1:N, "]")) %>% ggplot(aes(names, mean, group = 1)) + geom_line() + geom_point() +
  theme_imene(rotation_x = T)



#### GUILLAUME DATA ####
raw_data <- read.csv("Données_Guillaume.csv", sep = ";", col.names = c("Age_ka", "var"))
raw_data <- raw_data %>% mutate(sd = sqrt(var), depth = 1:17)
N = dim(raw_data)[1]

T1 <- min(raw_data$Age_ka)-3*min(raw_data$sd)
T2 <- max(raw_data$Age_ka)+3*max(raw_data$sd)
#visualisation 
p <- raw_data %>% uncount(1000) %>% mutate(value = rnorm(n(), Age_ka, sd)) %>% 
  ggplot(aes(x = value, y = depth, group = depth)) + 
  geom_density_ridges(scale = 1,fill= pallet[1], alpha = 0.4) + 
  geom_point(data = raw_data, aes(Age_ka, depth), color = pallet[2], inherit.aes = T, size = 3) + theme_imene() + 
  geom_segment(aes(x = T1, xend = T2, y = raw_data$depth[1], yend = raw_data$depth[1]), 
               inherit.aes = F, color = pallet[3])
p


##### Uniform Order #####
plain_order <- nimbleModel(ordre, list(N = N), 
                           data = list(M = raw_data$Age_ka, tau = raw_data$sd, T1 =T1, T2 = T2 ),
                           inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
cord <-  compileNimble(plain_order)
conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
mcmc <- buildMCMC(conf_ord)
cmcmc <-  compileNimble(mcmc)


##### Beta-Rademacher #####
rademacher_beta <- nimbleModel(inv_beta, list(N=N, p = 0.5), 
                                data = list(M = raw_data$Age_ka, tau = raw_data$sd, T1 = T1, T2 = T2),
                                inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                             b = rbeta(N, 1, N)),
                                dimensions = list(v = N))

cord_beta <- compileNimble(rademacher_beta)
conford_beta <- configureMCMC(rademacher_beta, monitors = c("mu", "Z", "b"), thin = 5)
mcmc_beta <- buildMCMC(conford_beta)
cmcmc_beta <- compileNimble(mcmc_beta)

##### BR(pi) ####
rademacher_betai <- nimbleModel(beta_pi, list(N=N, alpha = 0, beta = 1), 
                                data = list(M = raw_data$Age_ka, tau = raw_data$sd, T1 = T1, T2 = T2),
                                inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                             b = rbeta(N, 1, N), p = runif(N)),
                                dimensions = list(v = N))

cord_betai <- compileNimble(rademacher_betai)
conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "p", "b"), thin = 5)
mcmc_betai <- buildMCMC(conford_betai)
cmcmc_betai <- compileNimble(mcmc_betai)

##### BR(beta_pi) ####
rademacher_betapi <- nimbleModel(br_pi_beta, list(N=N, alpha = 4), 
                                data = list(M = raw_data$Age_ka, tau = raw_data$sd, T1 = T1, T2 = T2),
                                inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                             b = rbeta(N, 1, N), p = runif(N)),
                                dimensions = list(v = N))

cord_betapi <- compileNimble(rademacher_betapi)
conford_betapi <- configureMCMC(rademacher_betapi, monitors = c("mu", "Z", "p", "b"), thin = 5)
mcmc_betapi <- buildMCMC(conford_betapi)
cmcmc_betapi <- compileNimble(mcmc_betapi)

##### MCMC run ##### 
samp_beta <- runMCMC(cmcmc_beta, niter = 30000, nburnin = 22000, nchains = 3,
                      progressBar = TRUE, samplesAsCodaMCMC = TRUE)
samp_betai <- runMCMC(cmcmc_betai, niter = 30000, nburnin = 22000, nchains = 3,
                     progressBar = TRUE, samplesAsCodaMCMC = TRUE)
samp_betapi <- runMCMC(cmcmc_betapi, niter = 30000, nburnin = 22000, nchains = 3,
                      progressBar = TRUE, samplesAsCodaMCMC = TRUE)


samp_order <-  runMCMC(cmcmc, niter = 26000, nburnin = 18000, nchains = 3,
                       progressBar = TRUE, samplesAsCodaMCMC = TRUE)

plot(samp_betapi)
#Step
get_step_br(samp_betapi, "ZB_betapi_rawData.pdf", 4)

#### Aggregate chain 
sampBeta <- aggregat_samp(samp_beta)[, (2*N+1): (3*N)]
sampOrder <- aggregat_samp(samp_order)
sampBetai <- aggregat_samp(samp_betai)[, (2*N+1): (3*N)]
sampBetapi <- aggregat_samp(samp_betapi)[, (2*N+1): (3*N)]


sampZ <- aggregat_samp(samp_beta)[, 0:N]
apply(sampZ, 2,mean)

ic_o <- data.frame(t(apply(sampOrder, 2, interval_credible)[1:2, ]), 
                   ages = factor(colnames(sampOrder), levels = colnames(sampOrder)), 
                   mesures = raw_data$Age_ka, model = rep("order", N))

ic <- data.frame(t(apply(sampBeta, 2, interval_credible)[1:2, ]), 
                 ages = factor(colnames(sampBeta), levels = colnames(sampOrder)), 
                 mesures = raw_data$Age_ka, model = rep("BR", N))

ic_BRpi <- data.frame(t(apply(sampBetai, 2, interval_credible)[1:2, ]), 
                 ages = factor(colnames(sampBetai), levels = colnames(sampOrder)), 
                 mesures = raw_data$Age_ka, model = rep("BR(pi)", N))

ic_betapi <- data.frame(t(apply(sampBetapi, 2, interval_credible)[1:2, ]), 
                      ages = factor(colnames(sampBetapi), levels = colnames(sampOrder)), 
                      mesures = raw_data$Age_ka, model = rep("BR(betapi)", N))

ic %>% bind_rows(ic_o) %>% bind_rows(ic_BRpi) %>% bind_rows(ic_betapi) %>% 
  ggplot(aes(x = ages, ymin = X1, ymax = X2, group = model, colour = model)) + 
  geom_linerange(position = position_dodge(width = .2)) +
  theme_imene(rotation_x = T) + geom_point(aes(ages, mesures), inherit.aes = F) 

ggsave("ic_comp_rawData.png")





#### Garbage ####
t = seq(-1, 1, length.out = 1000)
l = sapply(2:15, function(x,n) x**(2-n)/(1-x), x= t)
l <- cbind(l, 2:15)
as.data.frame(l) %>% mutate(t = t) %>% pivot_longer(!V15) %>% ggplot(aes(x = t, y = value, group = V15))

plot(t, l[, 1], type= 'l', ylim = c(-2000, 2000))
for (i in 3:15) {
  lines(t, l[, i], col = i)
}

