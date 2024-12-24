############################################################@
source("Code_nimble.R") #Charge les packages Nimble rcarbon
library(tidyverse)
library(readxl) #pour importer AgeOSL depuis le format excel
library(ggridges) #pour les graphes des densités 
library(ArchaeoPhases)
source("utils.R")
############################################################@

#### Increasing the number of samples with assumptions that all the ages are equal : delta = 0 (looking only for the oldest ?)

n = 5#seq(5,100, 5)

### import prior properties 
path = "./BR_priors/sim_n10.csv"
prior = read.csv(path, row.names = 1)
prior_order = read.csv("./BR_priors/sim_order_n10.csv", row.names = 1)


###-----------------------------------@
delta = 0 
sigma = 10 #uncertainty 

graph_ic <- list()
graph_pi <- list()
median_pi <- list()
traj_pi <- list() #median
#presumely the oldest ages for both models  
ic_order <- c()
ic_br <- c()
T1 <- 1000-3*sigma 
T2 <- 1000+3*sigma

for (j in seq_along(n)) {
  sim <- simulated_data(n[j], delta, sigma, 123) #number of samples growing 
  
  
  #### Applying the MCMC algorithm for BR and order bayesian model ####
  
  #### Uniform Order ###
  plain_order <- nimbleModel(ordre, list(N = n[j]), data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ),
                             inits = list(e = rexp(n[j]+1)), dimensions = list(v =n[j], mu =n[j]))
  cord <-  compileNimble(plain_order)
  conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
  mcmc <- buildMCMC(conf_ord)
  cmcmc <-  compileNimble(mcmc)
  
  samp_order <-  runMCMC(cmcmc, niter = 20000, nburnin = 12000, nchains = 3,
                         progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
  
  #### Beta-Rademacher ###
  rademacher_betai <- nimbleModel(beta_pi, list(N=n[j], alpha = 0, beta = 1), 
                                  data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                                  inits = list(e = rexp(n[j]+1), Z1 = rbinom(n[j], 1, .5), 
                                               b = rbeta(n[j], 1, n[j]), p = runif(n[j])),
                                  dimensions = list(v = n[j]))
  
  cord_betai <- compileNimble(rademacher_betai)
  conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "p", "b"), thin = 5)
  mcmc_betai <- buildMCMC(conford_betai)
  cmcmc_betai <- compileNimble(mcmc_betai)
  
  
  samp_betai <- runMCMC(cmcmc_betai, niter = 20000, nburnin = 12000, nchains = 3,
                        progressBar = TRUE, samplesAsCodaMCMC = TRUE)

  
  
  ## aggregating chains 
  samp_pi <- aggregat_samp(samp_betai)[, (3*n[j] + 1): (4*n[j])] #p_1 ... P_n(j)
  sampBetai <- aggregat_samp(samp_betai)[, (2*n[j]+1):(3*n[j])] #ages from 21-30
  sampOrder <- aggregat_samp(samp_order)
  
  ic_o <- data.frame(t(apply(sampOrder, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampOrder), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("order", n[j]))
  
  ic <- data.frame(t(apply(sampBetai, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampBetai), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("BR", n[j]))
  
  ic_pi <- data.frame(t(apply(samp_pi, 2, interval_credible)[1:2, ]), ages = factor(colnames(sampBetai), levels = colnames(sampOrder)))
  
  graph_ic[[j]] <- ic %>% bind_rows(ic_o) %>% ggplot(aes(x = ages, ymin = X1, ymax = X2, group = model, colour = model)) + geom_linerange(position = position_dodge(width = .2)) +theme_imene() + geom_point(aes(ages, mesures), inherit.aes = F) +
    theme(axis.text.x = element_text(angle = 45))
  
  graph_pi[[j]] <- ic_pi %>% ggplot(aes(x = ages, ymin = X1, ymax = X2 )) + geom_linerange() +theme_imene(rotation_x = T) +
    labs(title = paste0("Ic p_i for ", n[j], " samples"))
  
  traj_pi[[j]] <- data.frame(Median = apply(samp_pi, 2, median), names = ic_pi$ages) %>% ggplot(aes(names, Median, group = 1)) + 
    geom_point(size = 1) +
    geom_line() + theme_imene(rotation_x = T)
  
  #obs of the shift property
  older_order <- sampOrder[, n[j] ] -1000 #centered 
  older_BR <- sampBetai[, n[j] ] - 1000
  
  ic_order <- rbind(ic_order, interval_credible(older_order)[1:2])
  ic_br <- rbind(ic_br, interval_credible(older_BR)[1:2])
  
} ##### END OF LOOP

#----------------------------------------------------------------------------------@
##### Creating the graphics for the latest age (comparaison with order) ########


ic_order <- data.frame(inf = ic_order[,1], sup = ic_order[,2], N = factor(n, levels = as.character(n)), model = rep("order", length(n)))
ic_br <- data.frame(inf = ic_br[,1], sup = ic_br[,2], N = factor(n, levels = as.character(n)), model = rep("BR", length(n)))
ic_max <- ic_order %>% bind_rows(ic_br)
head(ic_max)
ic_max %>% ggplot(aes(x = N,ymin = inf, ymax = sup, group = model, colour = model )) + geom_linerange(position = position_dodge(width = .3)) + ylab("IC(95%)") + xlab("number of samples N") + theme_imene() +
  theme(axis.text.x = element_text(angle = 45))

#### saving graphs 
path_graph <- "SimuSteierRom"
l = length(n)
for (j in 1:(l-1)) {
  print(graph_ic[[j]] )
  # ggsave(paste0(path_graph, "/Growing_n/n_", n[j], ".png"), width = 13, height = 7)
}

for (j in 1:(l-1)) {
  print(traj_pi[[j]] )
  # ggsave(paste0(path_graph, "/Growing_n/pi_n", n[j], ".png"), width = 13, height = 7)
}

traj_pi[[1]]


data.frame(Median = apply(samp_pi, 2, median), names = ic_pi$ages) %>% ggplot(aes(names, Median, group = 1)) + geom_point(size = 1) +
  geom_line() + theme_imene(rotation_x = T)

plot(density(samp_pi[, 1]))
ggsave


graph_pi[[1]]

mcmc1 <- coda::mcmc(samp_betai$chain1[, 1:100]*samp_betai$chain1[, 101:200])
mcmc2 <- coda::mcmc(samp_betai$chain2[, 1:100] *samp_betai$chain2[, 101:200])
mcmc3 <- coda::mcmc(samp_betai$chain3[, 1:100] *samp_betai$chain3[, 101:200])

list_mcmc <- mcmc.list(mcmc1, mcmc2, mcmc3)

plot(list_mcmc)

summary(list_mcmc)

barplot(table(samp_betai$chain1[, 4]))

hist(samp_betai$chain1[, 101])


#--------------------------------------------------------------------------------------------------------------------@
##test du BR avec 2 échantillon 
N = 18
sim <- simulated_data(N, delta, sigma, 123) #number of samples growing
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
plot(samp_betai)
sampBetai <- aggregat_samp(samp_betai)[, (2*N+1):(3*N)]
sampBetai


plain_order <- nimbleModel(ordre, list(N = N), data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ),
                           inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
cord <-  compileNimble(plain_order)
conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
mcmc <- buildMCMC(conf_ord)
cmcmc <-  compileNimble(mcmc)

samp_order <-  runMCMC(cmcmc, niter = 20000, nburnin = 12000, nchains = 3,
                       progressBar = TRUE, samplesAsCodaMCMC = TRUE)

sampOrder <- aggregat_samp(samp_order)



ic_o <- data.frame(t(apply(sampOrder, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampOrder), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("order", N))

ic <- data.frame(t(apply(sampBetai, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampBetai), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("BR", N))

ic %>% bind_rows(ic_o) %>% ggplot(aes(x = ages, ymin = X1, ymax = X2, group = model, colour = model)) + geom_linerange(position = position_dodge(width = .2)) +theme_imene() + geom_point(aes(ages, mesures), inherit.aes = F)





M = data.frame(X = sapply(c(0, 5,10,15,20,25,30), function(n,x) exp(-n*x) * sin(2*n*x), x = seq(0,5, length.out = 100)) , x = seq(0,20, length.out = 100))
M %>% pivot_longer(!x) %>% ggplot(aes(x, value, group = value , color = value )) + geom_line()

plot(M$x, M$X.2, type = 'l', xlab = "x", ylab = "fn")
for(i in 2:5){
  lines(M$x, M[, i], col = i)
}
legend("topright", col = 2:5, lty = 1, legend = paste0("n=", c(0, 5,10,15,20,25,30)))











