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
n.iter.2s = 1000
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
## !!! To do: (1) Fix survdat, (2) Fix data.jags.jm, (3) Fix model_JM2.bug

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
n.iter.jm = 1000
n.burnin.jm = 2000
seed.jm = 852

## Set arguments for JAGS function
data.jags.jm = list(Nobs = nrow(longdat), N = N, Nsurv = nrow(survdatnew), y = longdat$y, longid = longdat$id, 
                    survid = survdatnew$id, event = survdatnew$pevent, period = survdatnew$period, 
                    start = survdatnew$start, end = survdatnew$end, nk = nk)
inits.jags.jm = list(BP = rep(0, N), mu.BP = 0, sd.BP = 15, log.BPSD = rep(2, N), eta = rep(-5, nk), alpha1 = 0,
                    alpha2 = 0, mu.log.BPSD = 2, sd.log.BPSD = 0.5, rho=0, .RNG.seed=seed.jm)
model.jags.jm = "R/model_JM2.bug"
parameters.to.save.jm = c("eta", "alpha1", "alpha2", "mu.BP", "sd.BP", "mu.log.BPSD", "sd.log.BPSD")

# Fit joint model using JAGS  
jagsfit.jm2 <- jags(data = data.jags.jm, inits = list(inits.jags.jm), parameters.to.save = parameters.to.save.jm,
                   model.file = model.jags.jm, n.chains = 1, n.iter = n.iter.jm, n.burnin = n.burnin.jm, 
                   n.thin = 1, DIC = FALSE)
print(jagsfit.jm2)


  
  









	###########################
	## JOINT MODEL  ##
	
	nk <- 5     # JB: nk is the number of knots in the baseline hazard

# Select baseline hazard knots at time of death quantiles
k <- quantile(Time[Status==1],probs=seq(0,1,1/nk))
k[1] <- 0
k[(nk+1)] <- 20

# Transform the survival data so that each individual contributes one row for each time period (for a piecewise basline hazard)
survdat <- data.frame(id=1:N,Time,Status)
# RP: repeat each row 4 times:
survdat <- survdat[rep(1:N,rep(nk,N)),]
# RP: add period: 1, 2, 3, 4, 1, 2, 3, 4, etc.:
survdat$period <- rep(1:nk,N)
# RP: add start: these are the quantiles represented in k, repeating (e.g. 0, 3.8, 8.7, 14.6, 20, 0, 3.8, etc.):
survdat$start <- k[survdat$period]
# RP: returns whatever is the minimum bw stime and the quantile of the following period - i.e., for a given row...
	# $start and $end indicate the start and end time of that chunk of time.
survdat$end <- pmin(k[survdat$period+1],survdat$Time)
# RP: prune out any rows with start times greater than end times (I guess this happens as some stimes occur before the later quantiles?)
survdat<- survdat[survdat$end>survdat$start,]
# RP: regenerate event: event == 0 if $end in that row is less than the final event time for that person ($stime), else == old value for event;
	# ...so this equals 1 if the current row is the final event, and this is less than the admin censoring period (20):
survdat$Status <- ifelse(survdat$end<survdat$Time,0,survdat$Status)

# JB: Fit the joint model without correlation
Nobs=N*n
nrows = nrow(survdat)
data_3 = list("Nobs" = Nobs,
	              "N" = N,
	              "nrow" = nrows,
	              "y" = measurements,
	              "id" = rep(1:N, n),
	              "survid" = survdat$id,
	              "event" = survdat$Status,
	              "period" = survdat$period,
	              "start" = survdat$start,
	              "end" = survdat$end,
			  "nk" = nk
	              )
inits_3 = list(list(BP = rep(0, N),
	                    mu.BP = 0,
	                    sd.BP = 15,
	                    log.BPSD = rep(2, N),
	                    eta = rep(-5, nk),
	                    alpha1 = 0,
	                    alpha2 = 0,
	                    mu.log.BPSD = 2,
	                    sd.log.BPSD = 0.5,
				  .RNG.name="base::Marsaglia-Multicarry",
				 .RNG.seed=seed
				)
	               )
print(i)
out_jm_nocor <- jags(data=data_3, inits=inits_3, parameters.to.save=c("eta", "alpha1", "alpha2", "mu.BP", "sd.BP", "mu.log.BPSD", "sd.log.BPSD"),
            model.file="jointmodel_3_jags.bug", n.chains=1) 

# JB: Fit the joint model with correlation
	nrows = nrow(survdat)
	data_4 = list("Nobs" = Nobs,
	              "N" = N,
	              "nrow" = nrows,
	              "y" = measurements,
	              "id" = rep(1:N, n),
	              "survid" = survdat$id,
	              "event" = survdat$Status,
	              "period" = survdat$period,
	              "start" = survdat$start,
	              "end" = survdat$end,
			  "nk" = nk
	              )
	inits_4 = list(list(BP = rep(0, N),
	                    mu.BP = 0,
	                    sd.BP = 15,
	                    log.BPSD = rep(2, N),
	                    eta = rep(-5, nk),
	                    alpha1 = 0,
	                    alpha2 = 0,
	                    mu.log.BPSD = 2,
	                    sd.log.BPSD = 0.5,
	                    rho=0,
				 .RNG.name="base::Mersenne-Twister",
				 .RNG.seed=seed
	                    )
	               )
	print(i)
	out_jm_cor <- jags(data = data_4,
	            inits = inits_4,
	            parameters.to.save = c("eta",
	                                   "alpha1",
	                                   "alpha2",
	                                   "mu.BP",
	                                   "sd.BP",
	                                   "mu.log.BPSD",
	                                   "sd.log.BPSD"
	                                   ),
	            model.file = "jointmodel_4_jags.bug",
	            n.chains = 1)

	## RP added up to here
	###########################
		
	    result = list( "sim.id" = c(result$sim.id, paste("SIM",i,"_" , Sys.time(),sep="")),
                   "model.1" = rbind(result$model.1, c(coxfit1$coefficients, # I do not show the result of these model in the report
                                                        coxfit1.linear$coefficients,
                                                        coxfit1.longi.1$coefficients,
                                                        coxfit1.longi.2$coefficients)) ,
                   "model.2" = rbind(result$model.2, c(coxfit2$coefficients,
                                                        coxfit2.linear$coefficients,
                                                        coxfit2.longi.1$coefficients,
                                                        coxfit2.longi.2$coefficients)),
                   "model.3" = rbind(result$model.3,
                                     out_jm_nocor$BUGSoutput$summary[, "mean"]),
			"model.4" = rbind(result$model.4,
                                     out_jm_cor$BUGSoutput$summary[, "mean"]),
			 "model.1.SE" = rbind(result$model.1.SE, c(sqrt(diag(coxfit1$var)), # I do not show the result of these model in the report
                                                        sqrt(diag(coxfit1.linear$var)),
                                                        sqrt(diag(coxfit1.longi.1$var)),
                                                        sqrt(diag(coxfit1.longi.2$var)))),
			 "model.2.SE" = rbind(result$model.2.SE, c(sqrt(diag(coxfit2$var)), # I do not show the result of these model in the report
                                                        sqrt(diag(coxfit2.linear$var)),
                                                        sqrt(diag(coxfit2.longi.1$var)),
                                                        sqrt(diag(coxfit2.longi.2$var)))),
			 "model.3.SE" = rbind(result$model.3.SE,
			                      out_jm_nocor$BUGSoutput$summary[, "sd"]),
			  "model.4.SE" = rbind(result$model.4.SE,
			                      out_jm_cor$BUGSoutput$summary[, "sd"]),
           
                   "n.event" = c(result$n.event, sum(Status)),
                   "rhohat" = c(result$rhohat, r$rho),
                   "true.betaBP" = c(result$true.betaBP, beta.BP ),
			 "true.betaBPV" = c(result$true.betaBPV, beta.BPV),
			 "true.rho" = c(result$true.rho, rho),
                   "parameter.n" = c(result$parameter.n, n) )
    return(result)
}



# Number of replications
R = 10

start.time.A <- Sys.time()
##################
# JB: set the seed for R
set.seed(240117)
result = NULL
# JB: set the seed for JAGS
seed <- 1
for(i in 1:R){
  system.time(  	
	result <- run(seed=seed) 
	)
   seed <- seed+1
}
colnames(result$model.1.SE) <- c("BP", "BP.linear", "BP.longi.1", "BP.longi.2")
colnames(result$model.2.SE) <- c("BP", "BPV", "BP.linear", "BPV.linear", "BP.longi.1", "BPV.longi.1", "BP.longi.2", "BPV.longi.2")
write.csv(result, "../../../output/for_paper/for_revision/simulations_Table1_n4.csv", row.names = F)
end.time.A <- Sys.time()



### Calculate simulation summary statistics
getresults <- function(infile,out,truelogHR.BP=0.02, truelogHR.BPV=0.05){

    	result = read.csv(infile)#
	print(dim(result))

# Calculate the RMSE
	truevals <- rep(c(truelogHR.BP,truelogHR.BPV),4)
    	sqerror <- sweep(result[,6:13],2,truevals) 
    	sqerror <- as.data.frame(sqerror^2)
	
    
    result = list(    
		"est" = list("BP" = result$model.2.BP, "BPV" = result$model.2.BPV,
                               "BP.linear" = result$model.2.BP.linear, "BPV.linear" = result$model.2.BPV.linear,
                               "BP.longi.1" = result$model.2.BP.longi.1, "BPV.longi.1" = result$model.2.BPV.longi.1,
                               "BP.longi.2" =result$model.2.BP.longi.2, "BPV.longi.2" =result$model.2.BPV.longi.2),
		"se" = list("BP" = result$model.2.SE.BP, "BPV" = result$model.2.SE.BPV,
                               "BP.linear" = result$model.2.SE.BP.linear, "BPV.linear" = result$model.2.SE.BPV.linear,
                               "BP.longi.1" = result$model.2.SE.BP.longi.1, "BPV.longi.1" = result$model.2.SE.BPV.longi.1,
                               "BP.longi.2" =result$model.2.SE.BP.longi.2, "BPV.longi.2" =result$model.2.SE.BPV.longi.2),
		"sqerr" = list("BP" = sqerror$model.2.BP, "BPV" = sqerror$model.2.BPV,
                               "BP.linear" = sqerror$model.2.BP.linear, "BPV.linear" = sqerror$model.2.BPV.linear,
                               "BP.longi.1" = sqerror$model.2.BP.longi.1, "BPV.longi.1" = sqerror$model.2.BPV.longi.1,
                               "BP.longi.2" =sqerror$model.2.BP.longi.2, "BPV.longi.2" =sqerror$model.2.BPV.longi.2
				),
		"n.event" = result$n.event
		)
 
	i <- 1
	for(name in names(result$est)){
            	est = result$est[[name]]
			se = result$se[[name]]
            	sqerr = result$sqerr[[name]]
			uci = est+qnorm(0.975)*se
			lci = est-qnorm(0.975)*se
			
			meanest = mean(est)
			sdest = sd(est)
			rmse = sqrt(mean(sqerr))
			coverage = mean(truevals[i]>=lci & truevals[i]<=uci)

	           	out = c(out, paste(round(meanest,4), " (",round(sdest,4), ")", sep =""), paste(round(rmse,4)), paste(round(coverage, 4)) )
            	i <- i+1
        	}
	out
}
    




