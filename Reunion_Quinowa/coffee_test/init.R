require(devtools)
install_github("Maarten14C/coffee")
library(coffee) 


### 
set.seed(123)
?strat
sim.strat()

s = seq(0, 1, length.out = 1000)
plot(s, (1-s)* s**5, type = 'l')
plot(s,  s**(1-6)/(1-s), type = 'l')
