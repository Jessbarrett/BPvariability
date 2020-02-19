#############################################################################################
# 18/2/20 Jessica K Barrett                                                                 #
#                                                                                           #
# Example code to simulate a dataset and fit naive, two-stage and joint models              #
# to estimate the association between blood pressure variability and                        #
# cardiovascular disease                                                                    #                                                  
#                                                                                           #
# For details of the methods see Barrett et al, Statistics in Medicine 2019; 38: 1855-1868. #
#                                                                                           #
#############################################################################################

# Load packages
library(survival) 
library(R2jags) 


### Simulation parameters :
  N = 1500             # Number of individuals
  n = 4                # Number of repeated measurements per individual
  mu.BP =  120         # Mean of random effects distribution for usual BP (random intercept)
  sd.BP = 15           # SD of random effects distribution for usual BP (random intercept)   
  mu.log.BPSD = 2       # Mean of random effects distribution for log BP standard deviation (BP variability)
  sd.log.BPSD = 0.5     # SD of random effects distribution for log BP standard deviation (BP variability)
  rho = 0.5            # Correlation of usual BP with log BP standard deviation (BP variability)
  alphaBP = 0.02      # HR of usual BP with time to first CVD event
  alphaBPSD = 0.05     # HR of BP standard deviation with time to first CVD event
  

### Simulate a dataset
# Set the seed for data generation
  seed = 240117
  
## Generate random effects
  BP = rnorm(N, mu.BP, sd.BP)
# Normal distribution for BPSD conditional on BP
  BPSD = exp(rnorm(N, mu.log.BPSD + sd.log.BPSD/sd.BP*rho*(BP - mu.BP) , sqrt((1-rho**2))*sd.log.BPSD ))
  
## Generate longitudinal data
  y = rnorm(N*n, BP , BPSD)
# Standardise observations
  y <- y-mean(y)
# Longitudinal dataframe
  id <- rep(1:N,n)
  longdat <- data.frame(id,y)
  longdat <- longdat[order(longdat$id),]
  
## Generate survival data
# Fix gamma0 to get approx 20% events at T=20:
# CDF of Weibull distn with shape=2 is F(t) = 1-exp(-bt^2), where b=exp(gamma0+alphaBP*BP+alphaBPSD*BPSD)
# So F(20)=0.2 -> b=-log(0.8)/400 
# So to achieve 20% events at population mean, gamma0=log(-log(0.8))/400 -0.02*120 - 0.05*exp(2) 
  gamma0 = log(-log(0.8)/400)-0.02*120-0.05*exp(2)
  HR = exp(gamma0 + alphaBP*BP + alphaBPSD*BPSD ) # HR = Hazard Ratio
  time <- rweibull(N, scale = HR^(-1/2), shape=2)
# Censor at T=20
  event <- as.numeric(time<=20)
  time <- ifelse(time<20,time,20)
# Survival dataframe
  id <- 1:N
  survdat <- data.frame(id,time,event)

  
### Fit naive model
survdat$BP.naive = tapply(longdat$y,longdat$id,mean)  
survdat$BPSD.naive = tapply(longdat$y,longdat$id,sd)
coxfit.naive <- coxph(Surv(time,event)~BP.naive+BPSD.naive,data = survdat)
summary(coxfit.naive)


### Fit two-stage model
# Set number of iterations, burn-in and seed for JAGS
n.iter.2s = 2000
n.burnin.2s = 1000
seed = 741
  
# Set arguments for JAGS function
data.jags.2s = list(Nobs = nrow(longdat), N = N, y = longdat$y, subject = longdat$id)
inits.jags.2s = list(BP = rep(0,N),  logBPSD = rep(2,N), muBP = 0, sdBP = 15, mulogBPSD = 2, sdlogBPSD = 0.5,
                    rho=0, .RNG.seed=seed)
model.jags.2s = "R/model_LMM2.bug"
parameters.to.save.2s = c("BP", "BPSD", "rho")

# Fit linear mixed effects model using JAGS  
jagsfit.lmm2 = jags(data = data.jags.2s, inits = list(inits.jags.2s), parameters.to.save = parameters.to.save.2s, 
                    model.file = model.jags.2s, n.chains = 1, n.iter = n.iter.2s, n.burnin = n.burnin.2s, 
                    n.thin = 1, DIC = FALSE)

# Extract posterior means for BP and BP standard deviation from JAGS model fit 
survdat$BP.lmm2 = jagsfit.lmm2$BUGSoutput$mean$BP
survdat$BPSD.lmm2 = jagsfit.lmm2$BUGSoutput$mean$BPSD

# Fit survival model
coxfit.lmm2 <- coxph(Surv(time,event)~BP.lmm2+BPSD.lmm2,data = survdat)
summary(coxfit.lmm2)


### Fit joint model

nk <- 15     # nk is the number of knots in the baseline hazard

## Select baseline hazard knots at time of death quantiles
k <- quantile(survdat$time[survdat$event==1],probs=seq(0,1,1/nk))
k[1] <- 0
k[(nk+1)] <- 20

## Transform the survival data so that each individual contributes one row for each time period (for a piecewise baseline hazard)
survdatnew <- data.frame(id=1:N,survdat$time,survdat$event)
names(survdatnew) <- c("id","time","event")
# Replicate each row 4 times:
survdatnew <- survdatnew[rep(1:N,rep(nk,N)),]
# Add a time period
survdatnew$period <- rep(1:nk,N)
# Add the start time for each time period (i.e. the quantiles in k)
survdatnew$start <- k[survdatnew$period]
# Add the end time for each time period, i.e. the min of the survival time and the start time of the next time period
survdatnew$end <- pmin(k[survdatnew$period+1],survdatnew$time)
# Remove any rows with start times greater than end times 
survdatnew<- survdatnew[survdatnew$end>survdatnew$start,]
# Define pevent to be event indicator for an event taking place in that row's time period
survdatnew$pevent <- ifelse(survdatnew$end<survdatnew$time,0,survdatnew$event)

## Set number of iterations, burn-in and seed for JAGS
n.iter.jm = 3000
n.burnin.jm = 2000
seed.jm = 852

## Set arguments for JAGS function
data.jags.jm = list(Nobs = nrow(longdat), N = N, Nsurv = nrow(survdatnew), y = longdat$y, longid = longdat$id, 
                    survid = survdatnew$id, event = survdatnew$pevent, period = survdatnew$period, 
                    start = survdatnew$start, end = survdatnew$end, nk = nk)
inits.jags.jm = list(BP = rep(0, N), muBP = 0, sdBP = 15, logBPSD = rep(2, N), eta = rep(-5, nk), alpha1 = 0,
                    alpha2 = 0, mulogBPSD = 2, sdlogBPSD = 0.5, rho=0, .RNG.seed=seed.jm)
model.jags.jm = "R/model_JM2.bug"
parameters.to.save.jm = c("eta", "alpha1", "alpha2", "muBP", "sdBP", "mulogBPSD", "sdlogBPSD")

# Fit joint model using JAGS  
jagsfit.jm2 <- jags(data = data.jags.jm, inits = list(inits.jags.jm), parameters.to.save = parameters.to.save.jm,
                   model.file = model.jags.jm, n.chains = 1, n.iter = n.iter.jm, n.burnin = n.burnin.jm, 
                   n.thin = 1, DIC = FALSE)
# Display posterior means, posterior standard deviations and quantiles
print(jagsfit.jm2)





