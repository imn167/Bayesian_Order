############################################################@
source("Code_nimble.R") #Charge les packages Nimble rcarbon
library(tidyverse)
library(readxl) #pour importer AgeOSL depuis le format excel
library(ggridges) #pour les graphes des densités 
library(ArchaeoPhases)
source("utils.R")
############################################################@
##### Tableau de données avec une inversion #####

simulated_data <- function(n_ages, delta, sigma,inverted = F, n_inv = 2,  half = 1000){
  ### simulation des données 
  Y <- rnorm(n_ages, sd = sigma)
  mean <- 1000  + (1:n_ages) - (N+1) /2
  M <-  Y + mean
  depth <- seq(200, 400, length.out = n)
  
  ####@
  if (inverted){
  set.seed(123)  
  inverted <- sample(1:n_ages, 2) #inversion de 7 et 8 
  values <- dates[inverted]
  dates[inverted] <- values[n_inv:1]
  }
  ####@
  # dates[5] <- dates[5] + dt/2
  
  dd <- data.frame(Depth = depth, Age = mean, std = rep(sigma, n_ages), mesure = M,
                   Names = factor(paste0("mu[", 1:n_ages, "]"), levels = paste0("mu[", 1:n_ages, "]")))
  return(dd)
  
}

lag = 2
Data <- simulated_data(2,20, lag, 2)

ggplot(data = Data, mapping = aes(Age, Depth, color = Names)) +geom_point(size = 3) + theme_minimal() +
  labs(title = "Jeu de données simulé", subtitle = "inversion mu[7] et mu[8] et écart de l'ordre pour mu[8] ")

N <- length(Data$Age)

Data$Age <- rep(5, N)

############################################################@


model_A <- nimbleModel(order_A, list(N = N), 
                         data = list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22),
                       inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .7), p = runif(N)), 
                         dimensions = list(v =N))
cord_A <- compileNimble(model_A)
conford_A <- configureMCMC(model_A, monitors = c("mu", "Z", "p"), thin = 5)
mcmc_A <- buildMCMC(conford_A)
cmcmc_A <- compileNimble(mcmc_A)

samp_A <- runMCMC(cmcmc_A, niter = 20000, nburnin = 12000, nchains = 3,
                     progressBar = TRUE, samplesAsCodaMCMC = TRUE)

path0 = "Order_On_Ai/"
pdf(paste0(path0, "output_model.pdf"))
plot(samp_A)

sampA <- rbind(samp_A$chain1[, 11:20], samp_A$chain2[, 11:20], samp_A$chain3[, 11:20])
apply(sampA, 2, mean)
Data[c(2,4)]
Z_A <- rbind(samp_A$chain1[, 1:10], samp_A$chain2[, 1:10], samp_A$chain3[, 1:10])



# ggplot(Zp1_mean) + geom_line(aes(x = age, y = mean, color = p)) + geom_point(aes(x = age, y = mean))


mu_names = paste0("mu[", 1:10, "]")

dt_Z <- data.frame(Z_mean = apply(Z_A, 2, mean), Age = factor(mu_names, levels = mu_names, ordered = TRUE))

ggplot(dt_Z) + geom_line(aes(x = Age, y = Z_mean, group = 1)) + geom_point(aes(x = Age, y = Z_mean)) + theme_minimal()

dt <- as.data.frame(t(apply(sampA, 2, interval_credible)[-3,]))
dt <- dt %>% 
  mutate(Ages = factor(mu_names, levels = mu_names, ordered = TRUE), estimation = apply(sampA, 2, mean), Mean = Data$Age, 
         Depth = factor(Data$Depth, levels = as.character(Data$Depth)))

ggplot(dt, aes(x = Ages, ymin = V1, ymax = V2)) +
  geom_linerange(position = position_dodge(.2), color = "#0000FFA0") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + labs(title = "Credibility Interval at 95%")+
  geom_point(mapping = aes(x = Ages, y = Mean), data = dt,
             inherit.aes = FALSE, alpha = .85, size = 3, color = "green", fill = "white") +
  geom_point(mapping = aes(x =Ages, y = estimation), inherit.aes = F, color = "#0000FFA0", size = 3)


ttADependant <- datafr(sampA, Data,  Data$Depth, 10, 1, "estimate", "Depth")

### Tracé des densités 
Adependant_density <- ggplot(ttADependant, aes(x = estimate, y = Depth, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.16, 0.84), scale = 1, alpha =.3,
    rel_min_height = 1e-8)  +
  geom_point(data = dt, 
             aes(x = estimation, y = Depth), inherit.aes = FALSE, alpha = .85, size = 4, col = "#C733FF", shape = 7) +
  geom_point(data = dt, mapping = aes(x = Mean, y =Depth), inherit.aes = F, alpha = .85, size = 4, col = "green") +
  scale_fill_manual(
    name = "credibility interval at 68%", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0"),
    labels = c("(0, 0.16]","(0.16, 0.84)", "(0.84, 1]")
  ) +theme_ridges() +  theme(legend.position = c(.8, .2), legend.text.align = 0,
                             panel.grid = element_line(colour = "#A0DDAFA0"), 
                             legend.background = element_rect(fill = "#AADDDDA0")) + 
  labs(title = "Uniform order with inversion of epsilon between ages", 
       x = "Estimation (ka)",
       y = "Depth (cm)")
Adependant_density
ggsave(paste0(path0, "Posterior_Density.pdf"))




dev.off()

#############################@

rademacher_beta <- nimbleModel(inv_beta, list(N=N), 
                               data = list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22, p = .5),
                               inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), b = rbeta(N, 1, N)),
                               dimensions = list(v = N))

cord_beta <- compileNimble(rademacher_beta)
conford_beta <- configureMCMC(rademacher_beta, monitors = c("mu", "Z"), thin = 5)
mcmc_beta <- buildMCMC(conford_beta)
cmcmc_beta <- compileNimble(mcmc_beta)


cord_beta$setData(list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 20, p = .9))

samp_beta <- runMCMC(cmcmc_beta, niter = 20000, nburnin = 12000, nchains = 3,
                  progressBar = TRUE, samplesAsCodaMCMC = TRUE)
plot(samp_beta)




sampBeta <- rbind(samp_beta$chain1[, 11:20], samp_beta$chain2[, 11:20], samp_beta$chain3[, 11:20])

Z_beta <- rbind(samp_beta$chain1[, 1:10], samp_beta$chain2[, 1:10], samp_beta$chain3[, 1:10])

mu_names = paste0("mu[", 1:10, "]")
mu_names

dt_Z <- data.frame(Z_mean = apply(Z_beta, 2, mean), Age = factor(mu_names, levels = mu_names, ordered = TRUE))

ggplot(dt_Z) + geom_line(aes(x = Age, y = Z_mean, group = 1)) + geom_point(aes(x = Age, y = Z_mean)) + theme_minimal()

dt <- as.data.frame(t(apply(sampBeta, 2, interval_credible)[-3,]))
dt <- dt %>% 
  mutate(Ages = factor(mu_names, levels = mu_names, ordered = TRUE), estimation = apply(sampBeta, 2, mean), Mean = Data$Age, true = rep('tru', N))

ggplot(dt, aes(x = Ages, ymin = V1, ymax = V2)) +
  geom_linerange(position = position_dodge(.2), color = "#0000FFA0") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + labs(title = "Credibility Interval at 95%")+
  geom_point(mapping = aes(x = Ages, y = Mean), data = dt,
             inherit.aes = FALSE, alpha = .85, size = 3, color = "green", fill = "white") +
  geom_point(mapping = aes(x =Ages, y = estimation), inherit.aes = F, color = "#0000FFA0", size = 3)



#####################@

rademacher_betai <- nimbleModel(beta_pi, list(N=N, alpha = 0, beta = 1), 
                                data = list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22),
                               inits = list(e = rexp(N+1), Z1 = rbinom(N, 1, .5), 
                                            b = rbeta(N, 1, N), p = runif(N)),
                               dimensions = list(v = N))

cord_betai <- compileNimble(rademacher_betai)
conford_betai <- configureMCMC(rademacher_betai, monitors = c("mu", "Z", "p", "b"), thin = 5)
mcmc_betai <- buildMCMC(conford_betai)
cmcmc_betai <- compileNimble(mcmc_betai)

# cord_betai$setData(list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22))

samp_betai <- runMCMC(cmcmc_betai, niter = 20000, nburnin = 12000, nchains = 3,
                     progressBar = TRUE, samplesAsCodaMCMC = TRUE)

path0 = "Rademacher_pi/"

pdf(paste0(path0, "output_model.pdf"))
plot(samp_betai)



colnames(samp_betai$chain1)

sampb <- rbind(samp_betai$chain1[, 11:20], samp_betai$chain2[, 11:20], samp_betai$chain3[, 11:20])
apply(sampb, 2, mean)
sampBetai <- rbind(samp_betai$chain1[, 21:30], samp_betai$chain2[, 21:30], samp_betai$chain3[, 21:30])
apply(sampBetai, 2, mean)

Z_betai <- rbind(samp_betai$chain1[, 1:10], samp_betai$chain2[, 1:10], samp_betai$chain3[, 1:10])
p_betai <- rbind(samp_betai$chain1[, 31:40], samp_betai$chain2[, 31:40], samp_betai$chain3[, 31:40])


mu_names = paste0("mu[", 1:10, "]")

dt_Z <- data.frame(Z_mean = apply(Z_betai, 2, mean), Age = factor(mu_names, levels = mu_names, ordered = TRUE))
dt_p <- data.frame(p_mean = apply(p_betai, 2, mean), Age = factor(mu_names, levels = mu_names, ordered = TRUE))

ggplot(dt_p) + geom_line(aes(x = Age, y = p_mean, group = 1)) + geom_point(aes(x = Age, y = p_mean)) + theme_minimal()

ggplot(dt_Z) + geom_line(aes(x = Age, y = Z_mean, group = 1)) + geom_point(aes(x = Age, y = Z_mean)) + theme_minimal()

dt <- as.data.frame(t(apply(sampBetai, 2, interval_credible)[-3,]))
dt <- dt %>% 
  mutate(Ages = factor(mu_names, levels = mu_names, ordered = TRUE), estimation = apply(sampBetai, 2, mean), Mean = Data$Age, 
         Depth = factor(Data$Depth, levels = as.character(Data$Depth)))

ggplot(dt, aes(x = Ages, ymin = V1, ymax = V2)) +
  geom_linerange(position = position_dodge(.2), color = "#0000FFA0") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + labs(title = "Credibility Interval at 95%")+
  geom_point(mapping = aes(x = Ages, y = Mean), data = dt,
             inherit.aes = FALSE, alpha = .85, size = 3, color = "green", fill = "white") +
  geom_point(mapping = aes(x =Ages, y = estimation), inherit.aes = F, color = "#0000FFA0", size = 3)

ttRademacherpi <- datafr(sampBetai, Data,  Data$Depth, 10, 1, "estimate", "Depth")

### Tracé des densités 
Adependant_density <- ggplot(ttRademacherpi, aes(x = estimate, y = Depth, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.16, 0.84), scale = 1, alpha =.3,
    rel_min_height = 1e-8)  +
  geom_point(data = dt, 
             aes(x = estimation, y = Depth), inherit.aes = FALSE, alpha = .85, size = 4, col = "#C733FF", shape = 7) +
  geom_point(data = dt, mapping = aes(x = Mean, y =Depth), inherit.aes = F, alpha = .85, size = 4, col = "green") +
  scale_fill_manual(
    name = "credibility interval at 68%", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0"),
    labels = c("(0, 0.16]","(0.16, 0.84)", "(0.84, 1]")
  ) +theme_ridges() +  theme(legend.position = c(.8, .2), legend.text.align = 0,
                             panel.grid = element_line(colour = "#A0DDAFA0"), 
                             legend.background = element_rect(fill = "#AADDDDA0")) + 
  labs(title = "Ai = A(i-1) + Zi * Bi", 
       x = "Estimation (ka)",
       y = "Depth (cm)")
Adependant_density
ggsave(paste0(path0, "Posterior_Density.pdf"))

















dev.off()


#Illustrer le problème de l'amplitude pour Z rademacher c'est une limite quand les valeurs sont 1 ou -1

#######@

gaussian <-  nimbleModel(gaussien_noise, list(N=N), 
                         data = list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22),
            inits = list(e = rexp(N+1), Z = rnorm(N, 0, 1), b = rbeta(N, 1, N), sigma = rinvgamma(N, 1,1), 
                         alpha = runif(1, 0, 1), beta = runif(1, 0, 1)),
            dimensions = list(v = N))

cord_gaussian <- compileNimble(gaussian)
conf_gaussian <- configureMCMC(gaussian, monitors = c("mu", "Z", "sigma", "b", "alpha", "beta"), thin = 5)
mcmc_gaussian <-  buildMCMC(conf_gaussian)
cmcmc_gaussian <- compileNimble(mcmc_gaussian)

####@
# cord_gaussian$setData(list(M = Data$Age, tau = Data$std, T1 = 0, T2 = 22))
####@

samp_gaussian <- runMCMC(cmcmc_gaussian, niter = 50000, nburnin = 32000, nchains = 3,
                         progressBar = TRUE, samplesAsCodaMCMC = TRUE)

path0 = "GaussianNoise/"
pdf(paste0(path0, "output_model_unif.pdf"))
plot(samp_gaussian)



colnames(samp_gaussian$chain1)

sampGaussian <- rbind(samp_gaussian$chain1[, 23:32], samp_gaussian$chain2[, 23:32], 
                      samp_gaussian$chain3[, 23:32])
apply(sampGaussian, 2, mean)
Data$Age

Z_gaussian <- rbind(samp_gaussian$chain1[, 1:10], samp_gaussian$chain2[, 1:10], samp_gaussian$chain3[, 1:10])
apply(Z_gaussian, 2, mean)

dt_Z <- data.frame(Z_mean = apply(Z_gaussian, 2, mean), 
                   Age = factor(mu_names, levels = mu_names, ordered = TRUE))

ggz <- ggplot(dt_Z) + geom_line(aes(x = Age, y = Z_mean, group = 1)) + geom_point(aes(x = Age, y = Z_mean)) + theme_minimal() 
ggz
ggsave(paste0(path0, "Mean_Traj_Z.pdf"))

mu_names = paste0("mu[", 1:10, "]")
mu_names
dt <- as.data.frame(t(apply(sampGaussian, 2, interval_credible)[-3,]))
dt <- dt %>% 
  mutate(Ages = factor(mu_names, levels = mu_names, ordered = TRUE), 
         estimation = apply(sampGaussian, 2, mean), Mean = Data$Age, 
         Depth = factor(Data$Depth, levels = as.character(Data$Depth)))

gg_mu <- ggplot(dt, aes(x = Ages, ymin = V1, ymax = V2)) +
  geom_linerange(position = position_dodge(.2), color = "#0000FFA0") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) +
  geom_point(mapping = aes(x = Ages, y = Mean), data = dt,
             inherit.aes = FALSE, alpha = .85, size = 3, color = "green" ) +
  geom_point(mapping = aes(x =Ages, y = estimation), inherit.aes = F, color = "#0000FFA0", size = 3) + 
  labs(title = "Gaussian Noise model", subtitle = paste0("lag between order of ", lag))
gg_mu
ggsave(paste0(path0, "Age_estimation.pdf"))
Difference <-  data.frame(Age = mu_names, DiffProb = (Data$Age-apply(sampBetai, 2, mean)), 
                        GaussianNoise = (Data$Age - apply(sampGaussian, 2, mean)), 
                        A = (Data$Age - apply(sampA, 2, mean)))

gg_diff <- ggplot(data = Difference) + geom_line(aes(Age, DiffProb, group= 1), color = 'blue') +
  geom_point(mapping = aes(x = Age, y = DiffProb), color = "blue") +
  geom_line( aes(x = Age, y = GaussianNoise, group = 1), color = 'green') +
  geom_point(mapping = aes(Age, GaussianNoise), inherit.aes = F, color = "green") +
  geom_line( aes(x = Age, y = A, group = 1), color = 'red') +
  geom_point(mapping = aes(Age, A), inherit.aes = F, color = "red")+ 
  theme_minimal() +
  labs(title = "Quadratic Error", subtitle = "Difference between Rademacher(pi) et N(0, taui)")
gg_diff
ggsave(paste0(path0, "Squared Error.pdf"))



#### Representation graphique
### Calcul de l'estimation bayesienne 
Estimation <- data.frame(Depth = Data$Depth, Estimate = apply(sampGaussian, 2, mean), Real = Data$Age)
Estimation$Depth <- factor(Estimation$Depth, levels = as.character(Data$Depth))

ttGaussian <- datafr(sampGaussian, Data,  Data$Depth, 10, 1, "estimate", "Depth")

### Tracé des densités 
Gaussien_density <- ggplot(ttGaussian, aes(x = estimate, y = Depth, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.16, 0.84), scale = 1, alpha =.3,
    rel_min_height = 1e-8)  +
  geom_point(data = Estimation, 
             aes(x = Estimate, y = Depth), inherit.aes = FALSE, alpha = .85, size = 4, col = "#C733FF", shape = 7) +
  geom_point(data = Estimation, mapping = aes(x = Real, y =Depth), inherit.aes = F, alpha = .85, size = 4, col = "green") +
  scale_fill_manual(
    name = "credibility interval at 68%", values = c("#0000FFA0", "#FFF0AAAA", "#0000FFA0"),
    labels = c("(0, 0.16]","(0.16, 0.84)", "(0.84, 1]")
  ) +theme_ridges() +  theme(legend.position = c(.8, .2), legend.text.align = 0,
                             panel.grid = element_line(colour = "#A0DDAFA0"), 
                             legend.background = element_rect(fill = "#AADDDDA0")) + 
  labs(title = "Ai = ui + Zi *Bi", subtitle = "Zi N(0, tau(i)^2)", 
       x = "Estimation (ka)",
       y = "Depth (cm)")
Gaussien_density
ggsave(paste0(path0, "Posterior_Density.pdf"))


dev.off()






##### Importation des données de guillaume ####
mesures = read.csv("Données_Guillaume.csv", sep = ";") %>% mutate(sd = sqrt(variance))
colnames(mesures) <- c("Age_ka", "variance", "sd")
head(mesures)
N = dim(mesures)[1]

#lauching the BetaRademacher
rademacher_betai <- nimbleModel(beta_pi, list(N=N, alpha = 0, beta = 1), 
                                data = list(M = mesures$Age_ka, tau = mesures$sd, T1 = 0, T2 = 100),
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
sampBetai <- aggregat_samp(samp_betai)[, 35:51] 
dim(sampBetai)

IC_Betai = data.frame(inf = apply(sampBetai, 2, interval_credible)[1, ], 
                      sup = apply(sampBetai, 2, interval_credible)[2, ], 
                      estimation = apply(sampBetai, 2, median), #MAP meilleur, median voir les proprièté  ? 
                      mesure = mesures$Age_ka, names =factor( colnames(sampBetai), levels = colnames(sampBetai)))


IC_Betai %>% ggplot(aes(x = names, ymin = inf, ymax = sup)) + 
  geom_linerange(position = position_dodge(width = 0.2)) + 
  geom_point(aes(names, mesure, color = "mesure")) + 
  scale_colour_manual(
    values = c("mesure" = pallet[3])) + labs(color = "Legend") + xlab("Samples") + ylab("Ages (ka)") +
  theme_imene() + theme(axis.text.x = element_text(angle = 45))


#--------------------------------------------------@
plain_order <- nimbleModel(ordre, list(N = N), data = list(M = mesures$Age_ka, tau = mesures$sd, T1 =0, T2 = 100 ),
                           inits = list(e = rexp(N+1)), dimensions = list(v =N, mu =N))
cord <-  compileNimble(plain_order)
conf_ord <-  configureMCMC(plain_order, monitors = "mu", thin = 5)
mcmc <- buildMCMC(conf_ord)
cmcmc <-  compileNimble(mcmc)

samp_order <-  runMCMC(cmcmc, niter = 20000, nburnin = 12000, nchains = 3,
                       progressBar = TRUE, samplesAsCodaMCMC = TRUE)
sampOrder <- aggregat_samp(samp_order)

IC_order = data.frame(inf = apply(sampOrder, 2, interval_credible)[1, ], 
                      sup = apply(sampOrder, 2, interval_credible)[2, ], 
                      estimation = apply(sampOrder, 2, median), #MAP meilleur, median voir les proprièté  ? 
                      mesure = mesures$Age_ka, names = factor( colnames(sampOrder), levels = colnames(sampOrder)))

write_csv(IC_Betai, "Reunion_Quinowa/IC_BRi.csv")
write_csv(IC_order, "Reunion_Quinowa/IC_Order.csv")

IC = (IC_order %>% mutate( model = rep('order', N))) %>% bind_rows(IC_Betai %>% mutate( model = rep("BR", N))) 

(IC_order %>% mutate( model = rep('order', N))) %>% 
  bind_rows(IC_Betai %>% mutate( model = rep("BR", N))) %>% 
  ggplot(aes(x = names, ymin = inf, ymax = sup, color = model)) + 
  geom_linerange(position = position_dodge(width = 0.2)) + theme_imene() +
  geom_point(aes(x = names, y = mesure, color = "Obs")) +  labs(color = 'legend') + xlab("Samples") + ylab("Ages (ka)") +
  scale_colour_manual(
    values = c("Obs" = pallet[3],"order" = pallet[4], "BR" = pallet[2])) + theme(axis.text.x = element_text(angle = 45))


write_csv(IC_order, "Reunion_Quinowa/IC_Comparaison.csv")




