############################################################@
source("Code_nimble.R") #Charge les packages Nimble rcarbon
library(tidyverse)
library(readxl) #pour importer AgeOSL depuis le format excel
library(ggridges) #pour les graphes des densités 
library(ArchaeoPhases)
source("utils.R")
############################################################@

#### Increasing the number of samples with assumptions that all the ages are equal : delta = 0 (looking only for the oldest ?)

n = seq(5,100, 5)

### import prior properties 
path = "./BR_priors/sim_n10.csv"
prior = read.csv(path, row.names = 1)
prior_order = read.csv("./BR_priors/sim_order_n10.csv", row.names = 1)


###-----------------------------------@
#### Initialisation ####
delta = 0 
sigma = 100 #uncertainty 

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
  
  
  #### Bayesian Model ####
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

  #### MCMC PLOT ####
  path <- paste0("~/Desktop/Testing_Priors/Stratigraphic Priors/graphics/expB/mcmc_plot/N", n[j], ".pdf")
  bz <- get_step_br(samp_betai, path, 4)
  
  #### aggregating chains ####
  samp_pi <- aggregat_samp(samp_betai)[, (3*n[j] + 1): (4*n[j])] #p_1 ... P_n(j)
  sampBetai <- aggregat_samp(samp_betai)[, (2*n[j]+1):(3*n[j])] # mu_1 ... mu_n(j)
  sampOrder <- aggregat_samp(samp_order)
  
  ic_o <- data.frame(t(apply(sampOrder, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampOrder), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("order", n[j]))
  
  ic <- data.frame(t(apply(sampBetai, 2, interval_credible)[1:2, ])-1000, ages = factor(colnames(sampBetai), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("BR", n[j]))
  
  
  graph_ic[[j]] <- ic %>% bind_rows(ic_o) %>% ggplot(aes(x = ages, ymin = X1, ymax = X2, group = model, colour = model)) + geom_linerange(position = position_dodge(width = .2)) +theme_imene() + geom_point(aes(ages, mesures), inherit.aes = F) +
    theme(axis.text.x = element_text(angle = 45))
  
  
  traj_pi[[j]] <- data.frame(Median = apply(samp_pi, 2, median), names = factor(paste0("pi[", 1:n[j], "]"), levels = paste0("pi[", 1:n[j], "]"))) %>% 
    ggplot(aes(names, Median, group = 1)) + 
    geom_point(size = 1) +
    geom_line() + theme_imene(rotation_x = T)
  
  #obs of the shift property
  older_order <- sampOrder[, n[j] ] -1000 #centered 
  older_BR <- sampBetai[, n[j] ] - 1000
  
  ic_order <- rbind(ic_order, interval_credible(older_order)[1:2])
  ic_br <- rbind(ic_br, interval_credible(older_BR)[1:2])
  
} #### END OF LOOP####

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

traj_pi[[10]]


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
N = 25
sigma = 10
delta = 0
T1 <- 1000-3.5*sigma 
T2 <- 1000+3.5*sigma
sim <- simulated_data(N, delta, sigma, 13789) #number of samples growing

p <- sim %>% uncount(1000) %>% mutate(value = rnorm(n(), Age, std)) %>% ggplot(aes(x = value, y = Depth, group = Depth)) + 
  geom_density_ridges(scale = 1,fill= pallet[1], alpha = 0.4) + 
  geom_point(data = sim, aes(Mesure, Depth), color = pallet[2], inherit.aes = T, size = 3) + theme_imene()
p <- p + geom_segment(aes(x= T1, xend = T2, y = sim$Depth[1], yend = sim$Depth[1]), linewidth = 0.5, color = pallet[3])
p
rademacher_betai <- nimbleModel(beta_pi, list(N=N, alpha = 0, beta = 1), 
                                data = list(M = sim$Mesure, tau = sim$std, T1 = T1, T2 = T2),
                                inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                             b = rbeta(N, 1, N), p = runif(N)),
                                dimensions = list(v = N))

cord_betai <- compileNimble(rademacher_betai)
conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "p", "b", "v"), thin = 5)
mcmc_betai <- buildMCMC(conford_betai)
cmcmc_betai <- compileNimble(mcmc_betai)


samp_betai <- runMCMC(cmcmc_betai, niter = 20000, nburnin = 12000, nchains = 3,
                      progressBar = TRUE, samplesAsCodaMCMC = TRUE)
plot(samp_betai)

L = list()
for (i in 1:3) {
  L[[i]] <- coda::mcmc(samp_betai[[i]][, 0:N] *samp_betai[[i]][, (N+1):(2*N)])
  colnames(L[[i]]) <- paste0("ZB[", 1:N, "]")
}
mcmc_list = coda::mcmc.list(L)

path <- "~/Desktop/Testing_Priors/Stratigraphic Priors/graphics/expB/mcmc_plot/N25.pdf"
pdf(path)
plot(mcmc_list)
dev.off()
sampBetai <- (aggregat_samp(samp_betai)[, (2*N+1):(3*N)] - T1)/(T2-T1)
sampV <- aggregat_samp(samp_betai)[, (3*N+1):(4*N)]

plain_order <- nimbleModel(ordre, list(N = N), data = list(M = sim$Mesure, tau = sim$std, T1 =T1, T2 = T2 ),
                           inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
cord <-  compileNimble(plain_order)
conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
mcmc <- buildMCMC(conf_ord)
cmcmc <-  compileNimble(mcmc)

samp_order <-  runMCMC(cmcmc, niter = 20000, nburnin = 12000, nchains = 3,
                       progressBar = TRUE, samplesAsCodaMCMC = TRUE)

sampOrder <- (aggregat_samp(samp_order) -T1)/(T2-T1)



ic_o <- data.frame(t(apply(sampOrder, 2, interval_credible)[1:2, ]), ages = factor(colnames(sampOrder), levels = colnames(sampOrder)), mesures = sim$Mesure-1000, model = rep("order", N))

ic <- data.frame(t(apply(sampBetai, 2, interval_credible)[1:2, ]), ages = factor(colnames(sampBetai), levels = colnames(sampBetai)), mesures = sim$Mesure-1000, model = rep("BR", N))

ic %>% bind_rows(ic_o) %>% ggplot(aes(x = ages, ymin = X1, ymax = X2, group = model, colour = model)) + 
  geom_linerange(position = position_dodge(width = .2)) +theme_imene(rotation_x = T) +
  geom_hline(yintercept = 0.5, linetype = 2)
  # geom_point(aes(ages, mesures), inherit.aes = F) 


path <- "~/Desktop/Testing_Priors/Stratigraphic Priors/graphics/expB/N25.png"
ggsave(path)



## pi 
samp_pi <- aggregat_samp(samp_betai)[, (3*N+1):(4*N)]
median_p <- data.frame(median = apply(samp_pi, 2, median), 
           p = factor(paste0("p[",1:N, "]"), levels = paste0("p[",1:N, "]"))) %>% 
  ggplot(aes(x = p, y = median, group = 1)) +
  geom_line() + geom_point() + theme_imene(rotation_x = T)
median_p

path <- "~/Desktop/Testing_Priors/Stratigraphic Priors/graphics/expB/P25.png"
ggsave(path)
p
curve(dbeta(x, 2, 2), from = 0, to= 1)
