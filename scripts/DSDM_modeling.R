library(depmixS4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)

setwd("C:\\Users\\thchen7\\Downloads\\coding")

### dataset ###
df3 <- readRDS("Dnase-seq_data.rds")

# Modeling Starts
finalModel = c()
minBIC = .Machine$double.xmax
set.seed(11111)
BIClist = c()

ntimes_trimmed <- rep(99, 7464)
replicates = 10
NumOfStates =9
maxit = 500         
tol = 1e-08        
crit = "relative"  

for(j in 1:replicates){
  #pacf0.2
  mod2 <- depmix(list(MGW~1+MGWlag1+MGWlag2+MGWlag3,
                      EP~1+EPlag1,
                      ProT~1+ProTlag1,
                      Roll~1,
                      HelT~1+HelTlag1),
                 data=df3,
                 transition = ~1,
                 nstates=NumOfStates,
                 family=list(gaussian(),
                             gaussian(),gaussian(),gaussian(),gaussian()), 
                 ntimes=ntimes_trimmed)
  
  fmod2 <- fit(mod2,
               emcontrol = em.control(maxit, tol, crit),
               verbose=TRUE)
  summary(fmod2)
  
  if(BIC(fmod2)<minBIC){
    finalModel = fmod2
    minBIC = BIC(fmod2)
  }
  BIClist = c(BIClist, BIC(fmod2))
}
saveRDS(finalModel,"Dnase-seq_model.rds")

stateSequence = posterior(finalModel)$state
