###### Code pour le modèle ########
library(nimble, nimbleCarbon)
# library(rcarbon)
library(docstring)
library(zipfR)

######################## Cumsum dans nimble#################
cumsum_nimble  <- nimbleFunction(
    run = function(x = double(1)) {
        returnType(double(1))
        n  <-  length(x)
        v  <- rep(1, n)
        v[1]  <- x[1]
        for (i in 2:n) {
           v[i]  <- v[i - 1] + x[i]
        }
        return(v)
    }
)
################ nimble_sort #######
nimble_sort <- nimbleFunction(
  run = function(x = double(1)){
    
    l = length(x)
    lx = l-1
    madeChange <- TRUE
    while (madeChange) {
      madeChange <- FALSE
      for (i in 1:lx) {
        i1 <- i + 1
        if (x[i] > x[i1]) { #if equal no need to change 
          tmp <- x[i]
          x[i] <- x[i1]
          x[i1] <- tmp
          madeChange <- TRUE
        }
      }
    }
    return(x)
    returnType(double(1))
    
  }
)

sort_x <- 
  
################ distribution shrinkage #######
dshrink  <- nimbleFunction(
  run = function(x = double(0), s = double(0), log = integer(0, default = 1)) {
    returnType(double(0))
     logProb  <- log(s^2 / (s^2 + x)^2) #la proba en echelle log
    if (log) {
      return(logProb)
    }else {
      return(exp(logProb))
    }
    return(logProb)
  }
)
rshrink  <- nimbleFunction(
  run = function(n = integer(0), s = double(0)) {
    returnType(double(0))
    u  <- runif(1) #on traite le cas de dim = 1 seulement 
    sigma2  <- (1 - u)^(-1) - 1
    return(s**2 * sigma2)
  }
)
# Si modifications de la distribution register
#registerDistributions( list(
 # dshrink = list(
  #  BUGSdist = "dshrink(s)",
   # range = c(0, Inf))))

############################################@
###### MixUnif ########
mix_unif  <- nimbleCode({
  #modele
  for (i in 1:N) {
    M[i] ~ dnorm(mean = mu[i], sd = tau[i])
    # proba poue le melange
    Z[i] ~ dbinom(size = 1, prob = p)
    #2 loi du mélange : pi2
    mu2[i] ~ dunif(min = T1, max = T2) #on fixe periode dans cst
  }
  #creation de l'ordre pour les uniformes
  for (i in 1:(N + 1)) {
    e[i] ~ dexp(rate = 1)
  }
  es <-  sum(e[])
  #premiere loi du melange : pi1
  mu1[] <-  cumsum_nimble(e[1:N]) / es
  #Loi de melange : loi a priori pour notre modele
  for (i in 1:N) {
    mu[i]  <- Z[i] * (mu1[i] * (T2 - T1) + T1) + (1 - Z[i]) * mu2[i]
  }

})
mix_unifC14  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = sqrt(tau[i]^2 + sigma[i]^2))
  }
  ### prior
  for (i in 1:N) {
       # proba poue le melange
    Z[i] ~ dbinom(size = 1, prob = p)
    #2 loi du mélange : pi2
    mu2[i] ~ dunif(min = T1, max = T2) #on fixe periode dans cst
  }
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  #uniform order
  es <-  sum(e[])
  mu1[] <-  cumsum_nimble(e[1:N]) / es

  for (i in 1:N) {
    # calibrated Ages
     A[i]  <- Z[i] * (mu1[i] * (T2 - T1) + T1) + (1 - Z[i]) * mu2[i]

     #uncalibrated ages
     mu[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14BP[])
     sigma[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14err[])
  }
})

############################################@

#### Reduc Shift ####

unif_shift_order <-  nimbleCode({
  
  #likelyhood 
  for (i in 1:N) {
    M[i] ~ dnorm(mu[i], sd = tau[i])
  }
  #Order according to the uniform shift 
  s ~ dunif(min = 0, max = (T2-T1))
  mu[1] ~ dunif(min = T1, max = (T2-s))
  mu[N] <- s + mu[1]
  for (i in 1:(N-2)) {
    u[i] ~ dunif(mu[1], mu[N])
  }
  
  mu[2:(N-1)] <- nimble_sort(u[])
  
})




############################################@

##### Ordre ########

ordre <- nimbleCode({
  #modèle
for (i in 1:N){
    M[i] ~ dnorm(mu[i], sd = tau[i])
  }
  #l'echatillon d'exponentielles
for (i in 1:(N + 1)){
    e[i] ~ dexp(1)
    }
es <-   sum(e[])
  #loi beta Be(i, n-i+1) :  echantillon des uniformes ordonnées
  v[]  <- cumsum_nimble(e[1:N]) / es
  mu[]  <- v[] * (T2 - T1) + T1
})

orderC14  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = sqrt(tau[i]^2 + sigma[i]^2))
  }
  ### prior
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  #uniform order
    # calibrated Ages
  A[] <-  (cumsum_nimble(e[1:N]) / es) * (T2 - T1) + T1

  for (i in 1:N) {
     #uncalibrated ages
     mu[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14BP[])
     sigma[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14err[])
  }
})
####################################################@
#####################################################@
#### Uniforme ####

uniform <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i] ~ dnorm(mean = mu[i], sd = tau[i])
    #prior
    mu[i] ~ dunif(min = T1, max = T2)
  }
})

uniformC14  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = sqrt(tau[i]^2 + sigma[i]^2))
    #prior
    # calibrated Ages
    A[i] ~ dunif(min = T1, max = T2)
  }
  for (i in 1:N) {
     #uncalibrated ages
     mu[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14BP[])
     sigma[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14err[])
  }
})

####################################################@
#### Partial Order ####
ineg_large  <- nimbleCode({
  for (i in 1:N) {
     M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
     Z[i]  ~ dbinom(size = 1, prob = p)
  }
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  v[]  <- cumsum_nimble(e[1:N]) / sum(e[])
  mu[1]  <- v[1] * (T2 - T1) + T1
  for (i in 2:N) {
     mu[i]  <- Z[i] * (v[i] * (T2 - T1) + T1) +
     (1 - Z[i]) * (v[i - 1] * (T2 - T1) + T1)
  }
})
####################################################@
#### Inversion ####
inv_epsilon  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:(N - 1)) {
    Z1[i]  ~ dbinom(size = 1, prob = p)
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
  }
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  mu[1]  <- v[1] * (T2 - T1) + T1
  for (i in 2:N) {
     mu[i]  <- (v[i-1] + Z[i - 1] * (v[i] - v[i - 1])) * (T2 - T1) + T1
  }
})
###======================================================================@

### BR p fixe ####
inv_beta  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    Z1[i]  ~ dbinom(size = 1, prob = p)
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
     mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})

inv_betaC14  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = sqrt(tau[i]^2 + sigma[i]^2))
  }
  for (i in 1:N) {
    Z1[i]  ~ dbinom(size = 1, prob = p)
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
     e[i]  ~ dexp(rate = 1)
  }
  #uniform order
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])

  for (i in 1:N) {
    # calibrated Ages
     A[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
     #uncalibrated ages
     mu[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14BP[])
     sigma[i]  <- interpLin(z = abs(A[i]), x = calBP[], y = C14err[])
  }
})


######### Inversion  Be(1,n) et variable latente Uniforme sur (-1,1) #######

inv_uniform  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    Z[i]  ~ dunif(min= -1, max = 1)
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})

########## Deependant des ages Ai* #########

order_A <- nimbleCode(
  {
    #model
    for (i in 1:N) {
      M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
    }
    for (i in 1:N) {
      p[i] ~ dunif(min = 0, max = 1)
      Z1[i]  ~ dbinom(size = 1, prob = p[i])
      Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    }
    for (i in 1:(N + 1)) {
      e[i]  ~ dexp(rate = 1)
    }
    v[] <- cumsum_nimble(e[1:N]) / sum(e[])
    mu[1] <-  v[1] * (T2 - T1) + T1
    for (i in 2:N) {
      mu[i]  <- ( mu[i-1] + Z[i] * (v[i] - v[i-1]) * (T2 - T1) + T1 ) 
    }
  }
)



##### BR pi Unif ######


beta_pi <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    p[i] ~ dunif(min = alpha, max = beta)
    Z1[i]  ~ dbinom(size = 1, prob = p[i])
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})

##### BR pi beta ####

br_pi_beta <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    p[i] ~ dbeta(shape1 = alpha, shape2 =  alpha) #having a 
    Z1[i]  ~ dbinom(size = 1, prob = p[i])
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})


######### Modele Hierarchic bruit gaussien ######

gaussien_noise <-  nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }

  alpha ~ dunif(min = 0.1, max = 3)
  beta ~ dunif(min = 0.1, max = 3)
  for (i in 1:N) {
    # alpha[i] ~ dgamma(shape = 0.01, rate = 0.01)
    # beta[i] ~ dgamma(shape = 0.01, rate = 0.01)
    sigma[i] ~ dinvgamma(shape = alpha, scale = beta)
    Z[i]  ~ dnorm(mean = 0, sd = sigma[i])
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})




####### MODEL CDF gaussienne**** #####


inv_betai  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    x[i] ~ dunif(min = (M[i] - tau[i]), max = (M[i] + tau[i]))
    Z[i]  <- (Z1[i] - .5) * 2 #{-1, 1}
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})



#################### Bruit avec une beta et Z une uniforme sur (-a,a) ############

uniform_beta  <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i]  ~ dnorm(mean = mu[i], sd = tau[i])
  }
  for (i in 1:N) {
    Z[i]  ~ dunif(min = -a, max = a)
    b[i]  ~ dbeta(shape1 = 1, shape2 = N)
  }
  for (i in 1:(N + 1)) {
    e[i]  ~ dexp(rate = 1)
  }
  v[] <- cumsum_nimble(e[1:N]) / sum(e[])
  for (i in 1:N) {
    mu[i]  <- (v[i] + Z[i] * b[i]) * (T2 - T1) + T1
  }
})

##################### Chronomodel ################
chronomodel  <- nimbleCode({
  for (i in 1:N) {
    M[i] ~ dnorm(mean = mu[i], sd = tau[i])
    mu[i] ~ dnorm(mean = theta[i], sd = sigma[i])
    #prior de variance 
    sigma2[i] ~ dshrink(s = tau[i])
    sigma[i]  <- sqrt(sigma2[i])
  }
  for(i in 1:(N + 1)){
    e[i]   ~ dexp(rate = 1)
  }
  #order 
  theta[]  <- (cumsum_nimble(e[1:N]) / sum(e[])) * (T2-T1) + T1
})


#####################################################@
############################################@
#### semblant d'ordre ####

semblant_ordre <- nimbleCode({
  for (i in 1:N) {
    M[i] ~ dnorm(mu[i], sd = tau[i])
    #prior
    v[i] ~ dbeta(shape1 = i, shape2 = N - i + 1)
    mu[i] <- (v[i] * (T2 - T1) + T1)
    }
})

##### Melange ordre et semblant ordre avec les Be ########
melange_beta  <- nimbleCode({
  #modele
  for (i in 1:N) {
    M[i] ~ dnorm(mean = mu[i], sd = tau[i])
    # proba poue le melange
    Z[i] ~ dbinom(size = 1, prob = p)
    #2 loi du mélange : pi2
    mu2[i] ~ dbeta(shape1 = i, shape2 = N - i + 1) #def sur (0,1)
  }
  #creation de l'ordre pour les uniformes
  for (i in 1:(N + 1)) {
    e[i] ~ dexp(rate = 1)
  }
  es <-  sum(e[])
  #premiere loi du melange : pi1
  mu1[] <-  cumsum_nimble(e[1:N]) / es
  #Loi de melange : loi a priori pour notre modele
  for (i in 1:N) {
    mu[i]  <- Z[i] * (mu1[i] * (T2 - T1) + T1) +
    (1 - Z[i]) * (mu2[i] * (T2 - T1) + T1)
  }

})


########=====================================================================================================
#### Ecriture du data frame pour ggridge ####
## fonction qui reconstruit les données simulées
#pour le graph de densitées sur ggridges  
datafr  <- function(sampmcmc, data, fc_level, N, n_grp, nom_var1, nom_var2) { # nolint
  #'  {Cette fonction renvoi un tableau de données de 2 colonnes
  #'  \n une qui contient les estimations \n et une qui contient les profondeurs
  #' @param n_grp : le numéro de colonne de la variable qualitative
  #' @param fc_level : vecteur des valeurs de la variable qualitatives
  #' @param nom_var1_nom_var2: chaine de caractere des noms des deux variables
  #' @param N: la taille de l'echantillon
  #' @param n_grp: la position de la variable qui décrit \n l'échantillon (nom, position dans l'ordre) dans le data frame}

  #creation du tableau vide
  s  <- c()
  n  <- dim(sampmcmc)[1]
  for (i in 1:N) {
     grp  <- rep(data[i, n_grp], n)
     s <- rbind(s, cbind(sampmcmc[, i], grp))
  }
  s  <- as.data.frame(s)
  colnames(s) <- c(nom_var1, nom_var2)
  s[, 2]  <- as.character(s[, 2])
  s[, 2]  <- factor(s[, 2], levels = as.character(fc_level))
  return(s)
}

######### Fonction de simulation pour les Ages #######
sim_eps  <- function(n_ages, p, method = "methode1") {
  A  <- c()
  e  <- rexp(n_ages + 1)
  u  <- cumsum(e[1:n_ages]) / sum(e) #ordre
  Z  <- (rbinom(n_ages - 1, 1, p) - .5) * 2 #rademacher
  A  <- c(u[1], u[2:n_ages] + Z * u[2:n_ages])
  return(A)
}
####@
#Simulation de la loi a priori inv Beta(1,n)
sim_invbeta  <- function(n_ages, p, T1, T2) {
  A  <- c()
  e  <- rexp(n_ages + 1)
  u  <- cumsum(e[1:n_ages]) / sum(e) #ordre
  b  <-  rbeta(n_ages, 1, n_ages)
  Z  <- (rbinom(n_ages, 1, p) - .5) * 2 #rademacher
  A  <- (u + Z * b) * (T2-T1) + T1
  return(A)
}
######@
#Simulation de la loi a prior du melange uniforme
sim_melange  <- function(n_ages, p, T1, T2) {
  e  <- rexp(n_ages + 1)
  u  <- cumsum(e[1:n_ages]) / sum(e) #ordre
  b  <-  runif(n_ages)
  Z  <- rbinom(n_ages, 1, p) #rademacher
  A  <- (Z * u + (1 - Z) * b) * (T2 - T1) + T1
  return(A)
}
########@
#Simulation de la loi a priori pour ordre
sim_ordre  <- function(n_ages, T1, T2) {
  e  <- rexp(n_ages + 1)
  theta  <- (cumsum(e[1:n_ages]) / sum(e)) * (T2 - T1) + T1#ordre
  return(theta)
}
########@
#Simulation de la loi a priori pour ordre
sim_equal  <- function(n_ages, p,T1, T2) {
  e  <- rexp(n_ages + 1)
  ordre  <- (cumsum(e[1:n_ages]) / sum(e)) * (T2 - T1) + T1 #ordre
  Z  <- rbinom(n_ages-1, 1, p)
  A  <-  c(ordre[1], Z * ordre[2:n_ages] + (1 - Z) * ordre[1:(n_ages - 1)])
  return(A)
}
######### fonction de correlation pour l'estimation de echantillon i
corr_esp  <- function(i, j, n_samp, p) {
  if (i == j) {
     return(1)
  }else{
    return(i * (n_samp + 1 - j) / (n_samp + 1) **2 / (n_samp + 2))
  }
}

####################################@@@
##### Densité InvBeta ########
dens <- function(t, p, n, i){
  #t appartient à (-1, 2)
  cst = n * factorial(n)/ (factorial(i-1) * factorial(n-i)) #cst commune 
  k= 0:(n-1) #valeur de la k pour la somme 
  #expression de la densité sur (-1,0)
  tneg <- t[which(t<0, arr.ind = TRUE)]
  M10 <-  (1-p)*colSums(sapply(tneg, function(x,p) (1+x)^p, n-1-k)* sapply(tneg+1, Ibeta, k+i+1, n-i+1) *choose(n-1, k) * (-1)^k)
  #expression de la densite sur (0,1)
  t01 <- t[which(t<=1 & t>=0, arr.ind = TRUE)]
  M01 <- colSums(sapply(t01, function(x,p) (1-x)^p, n-1-k) * sapply(t01, Ibeta, k+i+1, n-i+1) * choose(n-1, k))*p + 
    colSums(sapply(t01, function(x,p) (x)^p, n-1-k) * beta(i+1, n+k-i+1) * choose(n-1, k)) *(1-p)
  #expression de la densite sur (1,2)
  t12 <- t[which(t>1, arr.ind = TRUE)]
  M12 <- p*colSums(sapply(t12, function(x,p) (1-x)^p, n-1-k) * (beta(k+1+i, n-i+1) - sapply(t12-1, Ibeta, k+i+1, n-i+1)) * choose(n-1, k))
  d <- cst *c(M10, M01, M12)
  tt = c(tneg, t01, t12)
  r = list(t = tt, density = d)
  return(r)
  
}

#deuxieme façon de generer la densite
dens1 <- function(len, p, n, i){
  #t appartient à (-1, 2)
  cst = n * factorial(n)/ (factorial(i-1) * factorial(n-i)) #cst commune 
  k= 0:(n-1) #valeur de la k pour la somme 
  #expression de la densité sur (-1,0)
  tneg <- seq(-1,0, length.out = len)
  M10 <-  (1-p)*colSums(sapply(tneg, function(x,p) (1+x)^p, n-1-k)* sapply(tneg+1, Ibeta, k+i+1, n-i+1) *choose(n-1, k) * (-1)^k)
  #expression de la densite sur (0,1)
  t01 <- seq(0,1, length.out = len)
  M01 <- colSums(sapply(t01, function(x,p) (1-x)^p, n-1-k) * sapply(t01, Ibeta, k+i+1, n-i+1) * choose(n-1, k))*p + 
    colSums(sapply(t01, function(x,p) (x)^p, n-1-k) * beta(i+1, n+k-i+1) * choose(n-1, k)) *(1-p)
  #expression de la densite sur (1,2)
  t12 <- seq(1,2, length.out = len)
  M12 <- p*colSums(sapply(t12, function(x,p) (1-x)^p, n-1-k) * (beta(k+1+i, n-i+1) - sapply(t12-1, Ibeta, k+i+1, n-i+1)) * choose(n-1, k))
  d <- cst *c(M10, M01, M12)
  tt = c(tneg, t01, t12)
  r = list(t = tt, density = d)
  return(r)
  
}


####################################################@
#### Fonction pour graphes lines ####
plotline  <- function(data, estimation, q1, q2, mean = TRUE) {
  q  <- apply(estimation, 2, quantile, probs = c(q1, q2))
  qinf  <- q[1, ]
  qsup  <- q[2, ]
  if (mean) {
     moy  <- apply(estimation, 2, mean)
     data  <- dplyr::mutate(data, estimate = moy, qinf = qinf, qsup = qsup)
  } else {
      data  <- dplyr::mutate(data, qinf = qinf, qsup = qsup)
  }
  return(data)
}

#######@
##### avec contrainte loi normale #########

code4 <- nimbleCode({
  #modèle
  for (i in 1:N) {
    M[i] ~ dnorm(mu[i], sd = tau[i])
  }
  #Loi priori
  for (i in 1:N) {
     e[i] ~ dexp(1)
     }
  #echantillon ordonné suivant la loi gaussienne
  mu[]  <-  qnorm(1 - exp(-cumsum_nimble(e[] / (N + 1 - (1:N)))),
                  mean = moy, sd = sigmaCurve)
})
####################################################@
#### Modele sans contraite d'ordre ####

code6 <- nimbleCode({
  #model
  for (i in 1:N) {
    M[i] ~ dnorm(mu[i], sd = tau[i])
    #prior
    mu[i] ~ dnorm(mean = moy, sd = sigmaCurve)
  }
})
####################################################@
#### Modele de melange Gaussian ####
code5  <- nimbleCode({
  #model
  for (i in 1:N) {
     M[i] ~ dnorm(mean = mu[i], sd = tau[i])
     #prior 1 : sans ordre
     mu2[i] ~ dnorm(mean = moy, sd = sigmaCurve)
     # proba de melange
     Z[i] ~ dbinom(size = 1, prob = p)
     #exp pour ordre
     e[i] ~ dexp(1)
  }
  #prior 2 : ordre
  mu1[]  <-  qnorm(1 - exp(-cumsum_nimble(e[] / (N + 1 - (1:N)))),
                mean = moy, sd = sigmaCurve)
  #melange
  for (i in 1:N) {
     mu[i]  <- Z[i] * mu1[i] + (1 - Z[i]) * mu2[i]
  }

})