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

  

  
  
  n.iter = 1500
  n.burnin = 500
  

###  Functions to estimate BP variability from repeated measurements
# Naive model
linear = function(measurements, n){
  r = matrix(measurements, ncol = n)
  return(list("BP" = apply(r, 1, mean), "BPV" = apply(r, 1, sd)))
}

# Longitudinal model
longi = function(measurements, model.number, n, seed){ 
  # model.number=1 longitudinal model with no correlation between BPSD and BP
  # model.number=2 longitudinal model with  correlation between BPSD and BP
  # JB: seed is the seed for JAGS
  data.bugs = list(Nobs = N*n, N = N, y = measurements, subject = rep(1:N, n))
  inits.bugs = list(BP = rep(0,N),  mu.BP = 0, sd.BP = 15, log.BPV = rep(2,N), mu.log.BPV = 2, sd.log.BPV = 0.5)
  model.bugs = paste("model_", model.number ,".bug", sep = "")
  parameters.to.save = c("BP", "BPV")
  if(model.number == 2){inits.bugs$rho = 0; parameters.to.save = c("BP", "BPV", "rho")} # In the second model rho is taking into account.

  inits.bugs$.RNG.seed=seed
  if(model.number==1){
	inits.bugs$.RNG.name="base::Super-Duper"
	} else {
	inits.bugs$.RNG.name="base::Wichmann-Hill"
	}
  
  r = jags(
    data = data.bugs,
    inits = list(inits.bugs),
    parameters.to.save = parameters.to.save,
    model.file = model.bugs,
    # model.file = "jointmodel_0a.bug",
    n.chains = 1,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = 1,
    DIC = FALSE,
  )
  return(list("BP" = unname(r$BUGSoutput$summary[1:N,1]),
              "BPV" = unname(r$BUGSoutput$summary[(N+1):(2*N),1]),
              "rho" = unname(r$BUGSoutput$summary[2*N+(model.number-1),1])
              )
         )
}

### Function to simulate and analyse a dataset
run = function(beta.BP = 0.02, beta.BPV = 0.05, rho = 0.5, n = 4, seed){
    cat("n = ", n,", beta.BP = ", beta.BP,", beta.BPV = ", beta.BPV,", rho = ", rho,", i = ", i,"\n") 
    
    

    ### Longitudinal analysis
      ## stage 1 : estimations of BP and BPV based on mesurements - linear, longi.1, longi.2
    
    # RP: deriving 'naive' estimates of usual BP and BPV:
    r = linear(measurements, n=n)
      COX$BP.linear = r$BP
      COX$BPV.linear = r$BPV
    
    # RP: deriving esimates of usual BP and BPV based on LMM1 (correlation bw random effects fixed at 0):
    r = longi(measurements, 1, n=n, seed=seed)
      COX$BP.longi.1 = r$BP
      COX$BPV.longi.1 = r$BPV

    # RP: deriving esimates of usual BP and BPV based on LMM2 (allows for non-zero correlation bw random effects):    
    r = longi(measurements, 2, n=n, seed=seed)
      COX$BP.longi.2 = r$BP
      COX$BPV.longi.2 = r$BPV
 
      ## stage 2 : Cox regression - and recording it in result list
	coxfit1 <- coxph(surv ~ BP , data = COX)
	coxfit1.linear <- coxph(surv ~ BP.linear , data = COX) 
	coxfit1.longi.1 <- coxph(surv ~ BP.longi.1 , data = COX)
	coxfit1.longi.2 <- coxph(surv ~ BP.longi.2 , data = COX)
	coxfit2 <- coxph(surv ~ BP + BPV , data = COX)	
	coxfit2.linear <- coxph(surv ~ BP.linear + BPV.linear , data = COX)
	coxfit2.longi.1 <- coxph(surv ~ BP.longi.1 + BPV.longi.1 , data = COX)
	coxfit2.longi.2 <- coxph(surv ~ BP.longi.2 + BPV.longi.2 , data = COX)
	
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
    




