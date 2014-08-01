#' 
#' @title Runs a full ESPRESSO analysis
#' @description This function calls the functions required to run a full ESPRESSO analysis 
#'  where the model consists of an outcome (binary or continuous) determined by two interacting
#'  covariates (a SNP  and an environmental exposure)
#' @param simulation.params general parameters for the scenario(s) to analyse
#' @param pheno.params paramaters for the outcome variables
#' @param env1.params parameters for the 1st environmental determinant
#' @param env2.params parameters for the 2nd environmental determinant
#' @param scenarios2run the indices of the scenarios one wish to analyse
#' @return a summary table that contains both the input parameters and 
#' the results of the analysis
#' @export
#' @author Gaye A.
#' @examples {
#'   
#' # load the table that hold the input parameters; each of the table
#' # hold parameters for 4 scenarios:
#' # scenario 1: a binary outcome determined by two binary environmental exposures 
#' # and an interaction between the two
#' # scenario 2: a binary outcome determined by two continuous environmental exposures 
#' # and an interaction between the two
#' # scenario 3: a quantitative outcome determined by two binary environmental exposures 
#' # and an interaction between the two
#' # scenario 4: a quantitative outcome determined by two continuous environmental exposures 
#' # and an interaction between the two
#' data(simulation.params) 
#' data(pheno.params)
#' data(env1.params)
#' data(env2.params)
#' 
#' # run the function for the first two scenarios, two binomial models
#' run.espresso.ExE(simulation.params, pheno.params, env1.params, env2.params, scenarios2run=c(1,2))
#'
#' # run the function for the last two scenarios, two gaussian models
#' run.espresso.ExE(simulation.params, pheno.params, env1.params, env2.params, scenarios2run=c(3,4))
#' 
#' }
#'
run.espresso.ExE <- function(simulation.params=NULL, pheno.params=NULL, env1.params=NULL, 
                             env2.params=NULL, scenarios2run=1){

# IF AN INPUT FILE IS NOT SUPPLIED LOAD THE DEFAULT TABLES WARNING
if(is.null(simulation.params)){
  cat("\n WARNING!\n")
  cat(" No simulation parameters supplied\n")
  cat(" The default simulation parameters will be used\n")
  simulation.params <- data(simulation.params)
}

if(is.null(pheno.params)){
  cat("\n WARNING!\n")
  cat(" No outcome parameters supplied\n")
  cat(" The default outcome parameters will be used\n")
  pheno.params <- data(pheno.params)
}

if(is.null(env1.params)){
  cat("\n WARNING!\n")
  cat(" No environmental parameters supplied for the 1st exposure\n")
  cat(" The default environmental parameters will be used\n")
  env1.params <- data(env1.params)
}

if(is.null(env2.params)){
  cat("\n WARNING!\n")
  cat(" No environmental parameters supplied for the 2nd exposure\n")
  cat(" The default environmental parameters will be used\n")
  env2.params <- data(env2.params)
}

# MERGE INPUT FILES TO MAKE ONE TABLE OF PARAMETERS
s.temp1 <- merge(simulation.params, pheno.params)
s.temp2 <- merge(s.temp1, env1.params)
s.parameters <- merge(s.temp2, env2.params)


#----------LOAD SET UP UP INITIAL PARAMETERS------------#

# PRINT TRACER CODE EVERY Nth ITERATION
# THIS ENSURES THAT YOU CAN SEE IF THE PROGRAM GRINDS TO A HALT FOR SOME REASON (IT SHOULDN'T)
trace.interval <- 10


# CREATE UP TO 20M SUBJECTS IN BLOCKS OF 20K UNTIL REQUIRED NUMBER OF
# CASES AND CONTROLS IS ACHIEVED. IN GENERAL THE ONLY PROBLEM IN ACHIEVING THE
# REQUIRED NUMBER OF CASES WILL OCCUR IF THE DISEASE PREVALENCE IS VERY LOW
allowed.sample.size <- 20000000
block.size <- 20000


# DECLARE MATRIX THAT STORE THE RESULTS FOR EACH SCENARIO (ONE PER SCENARIO PER ROW)
output.file <- "output.csv"
output.matrix <- matrix(numeric(0), ncol=48)
column.names <- c(colnames(s.parameters), "exceeded.sample.size?","numcases.required", "numcontrols.required", 
                 "numsubjects.required", "empirical.power", "modelled.power","estimated.OR", "estimated.effect")
write(t(column.names),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
                 
#-----------LOOP THROUGH THE SCENARIOS - DEALS WITH ONE SCENARIO AT A TIME-------------

for(j in c(scenarios2run))
{

   # RANDOM NUMBER GENERATOR STARTS WITH SEED SET AS SPECIFIED 
   set.seed(s.parameters$seed.val[j])

   # SIMULATION PARAMETERS
   scenario.id <- s.parameters$scenario.id[j]         
   seed.val <- s.parameters$seed.val[j]               
   numsims <- s.parameters$numsims[j]                 
   numcases <- s.parameters$numcases[j]               
   numcontrols <- s.parameters$numcontrols[j]  
   numsubjects <- s.parameters$numsubjects[j]
   int.OR <- s.parameters$interaction.OR[j]
   int.efkt <- s.parameters$interaction.efkt[j]            
   baseline.OR <- s.parameters$RR.5.95[j]                
   pval <- s.parameters$p.val[j]                      
   power <- s.parameters$power[j]
   
   # OUTCOME PARAMETERS
   pheno.model <- s.parameters$pheno.model[j]
   pheno.mean <- s.parameters$pheno.mean[j]
   pheno.sd <- s.parameters$pheno.sd[j]
   disease.prev <- s.parameters$disease.prev[j]
   pheno.error <- c(1-s.parameters$pheno.sensitivity[j],1-s.parameters$pheno.specificity[j])
   pheno.reliability <- s.parameters$pheno.reliability[j]    

   
   # 1st ENVIRONMENTAL DETERMINANT PARAMETERS
   env1.model<- s.parameters$env1.model[j]
   env1.prev <-  s.parameters$env1.prevalence[j]    
   env1.OR <- s.parameters$env1.OR[j]
   env1.efkt <- s.parameters$env1.efkt[j]
   env1.mean <- s.parameters$env1.mean[j]
   env1.sd <- s.parameters$env1.sd[j]
   env1.low.lim <- s.parameters$env1.low.lim[j]
   env1.up.lim <- s.parameters$env1.up.lim[j]
   env1.error <- c(1-s.parameters$env1.sensitivity[j],1-s.parameters$env1.specificity[j])
   env1.reliability <- s.parameters$env1.reliability[j]

   # 2nd ENVIRONMENTAL DETERMINANT PARAMETERS
   env2.model<- s.parameters$env2.model[j]
   env2.prev <-  s.parameters$env2.prevalence[j]    
   env2.OR <- s.parameters$env2.OR[j]
   env2.efkt <- s.parameters$env2.efkt[j]
   env2.mean <- s.parameters$env2.mean[j]
   env2.sd <- s.parameters$env2.sd[j]
   env2.low.lim <- s.parameters$env2.low.lim[j]
   env2.up.lim <- s.parameters$env2.up.lim[j]
   env2.error <- c(1-s.parameters$env2.sensitivity[j],1-s.parameters$env2.specificity[j])
   env2.reliability <- s.parameters$env2.reliability[j]   
      

   # VECTORS TO HOLD BETA, SE AND Z VALUES AFTER EACH RUN OF THE SIMULATION
   beta.values <- rep(NA,numsims)
   se.values <- rep(NA,numsims)
   z.values<-rep(NA,numsims)


   # TRACER TO DETECT EXCEEDING MAX ALLOWABLE SAMPLE SIZE
   sample.size.excess <- 0

   # GENERATE AND ANALYSE DATASETS ONE AT A TIME 
   for(s in 1:numsims)            # s from 1 to total number of simulations
   {

      #----------------------------------GENERATE "TRUE" DATA-----------------------------#
      

      if(pheno.model == 0){ # UNDER BINARY OUTCOME MODEL
        # GENERATE CASES AND CONTROLS UNTILL THE REQUIRED NUMBER OF CASES, CONTROLS IS ACHIEVED 
        sim.data <- sim.CC.data.ExE(n=block.size, ncases=numcases, ncontrols=numcontrols, 
                                    max.sample.size=allowed.sample.size, pheno.prev=disease.prev,
                                    e1.model=env1.model, e1.prev=env1.prev, e1.mean=env1.mean, 
                                    e1.sd=env1.sd, e1.low.lim=env1.low.lim, e1.up.lim=env1.up.lim, 
                                    e1.OR=env1.OR, e2.model=env2.model, e2.prev=env2.prev, 
                                    e2.mean=env2.mean, e2.sd=env2.sd, e2.low.lim=env2.low.lim, 
                                    e2.up.lim=env2.up.lim, e2.OR=env2.OR, i.OR=int.OR, 
                                    b.OR=baseline.OR, ph.error=pheno.error)

        true.data <- sim.data$data
        

      }else{ # UNDER QUANTITATIVE OUTCOME MODEL
        # GENERATE THE SPECIFIED NUMBER OF SUBJECTS
        true.data <- sim.QTL.data.ExE(n=numsubjects,ph.mean=pheno.mean,ph.sd=pheno.sd,
                                      e1.model=env1.model, e1.efkt=env1.efkt, e1.prev=env1.prev, 
                                      e1.mean=env1.mean, e1.sd=env1.sd, e1.low.lim=env1.low.lim, 
                                      e1.up.lim=env1.up.lim,e2.model=env2.model, e2.efkt=env2.efkt, 
                                      e2.prev=env2.prev, e2.mean=env2.mean, e2.sd=env2.sd, 
                                      e2.low.lim=env2.low.lim, e2.up.lim=env2.up.lim,
                                      i.efkt=int.efkt, pheno.rel=pheno.reliability)
      }

      #------------SIMULATE ERRORS AND ADD THEM TO THE TRUE COVARIATES DATA TO OBTAIN OBSERVED COVARIATES DATA-----------#

      # ADD APPROPRIATE ERRORS TO PRODUCE OBSERVED GENOTYPES 
      observed.data <- get.observed.data.ExE(data=true.data, e1.error=env1.error,e1.model=env1.model,
                                             e1.prev=env1.prev,e1.sd=env1.sd,e2.error=env2.error,
                                             e2.model=env2.model,e2.prev=env2.prev,e2.sd=env2.sd,
                                             e1.reliability=env1.reliability,e2.reliability=env2.reliability)

      #--------------------------DATA ANALYSIS ----------------------------#

      glm.estimates <- glm.analysis.ExE(pheno.model, observed.data)
     
      beta.values[s] <- glm.estimates[[1]]
      se.values[s] <- glm.estimates[[2]]
      z.values[s] <- glm.estimates[[3]]
      
      # PRINT TRACER AFTER EVERY Nth DATASET CREATED
      if(s %% trace.interval ==0)cat("\n",s,"of",numsims,"runs completed in scenario",scenario.id)

   }
   cat("\n\n")

   #------------------------ SUMMARISE RESULTS ACROSS ALL SIMULATIONS---------------------------#

   # SUMMARISE PRIMARY PARAMETER ESTIMATES
   # COEFFICIENTS ON LOG-ODDS SCALE
   mean.beta <- round(mean(beta.values, na.rm=T),3)
   mean.se <- round(sqrt(mean(se.values^2, na.rm=T)),3)
   mean.model.z <- mean.beta/mean.se
   
 
   #---------------------------POWER AND SAMPLE SIZE CALCULATIONS----------------------#

   # CALCULATE THE SAMPLE SIZE REQUIRED UNDER EACH MODEL
   sample.sizes.required <- samplsize.calc(numcases, numcontrols, numsubjects, pheno.model, pval, power, mean.model.z)

   # CALCULATE EMPIRICAL POWER AND THE MODELLED POWER 
   # THE EMPIRICAL POWER IS SIMPLY THE PROPORTION OF SIMULATIONS IN WHICH
   # THE Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC
   # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
   power <- power.calc(pval, z.values, mean.model.z)


   #------------------MAKE FINAL A TABLE THAT HOLDS BOTH INPUT PARAMETERS AND OUTPUT RESULTS---------------#

   critical.res <- get.critical.results.ExE(j,pheno.model,env1.model,env2.model,sample.sizes.required,power$empirical,
                                        power$modelled,mean.beta)

   #  WHEN OUTCOME IS BINARY INFORM IF RECORD EXCEEDED MAXIMUM SAMPLE SIZE
   if(pheno.model==0){
     sample.size.excess <- sim.data$allowed.sample.size.exceeded
     if(sample.size.excess==1)
     {
       excess <- "yes"
       cat("\nTO GENERATE THE NUMBER OF CASES SPECIFIED AT OUTSET\n")
       cat("THE SIMULATION EXCEEDED THE MAXIMUM POPULATION SIZE OF ", allowed.sample.size,"\n")
     }else{
       excess <- "no"
     }
   }
   
   inparams <- s.parameters[j,]
   if(pheno.model==0){
      mod <- "binary"
      if(env1.model==0){
        inparams [c(6,8,14,15,18,22:26,29)] <- "NA"
        inputs <- inparams
      }else{
        if(env1.model==1){
          inparams [c(6,8,14,15,18,20,21,25:28)] <- "NA"
          inputs <- inparams
        }else{
          inparams [c(6,8,14,15,18,20,21,23,24,27,28)] <- "NA"
          inputs <- inparams          
        }
      }
      if(env2.model==0){
        inparams [c(6,8,14,15,18,33:37,40)] <- "NA"
        inputs <- inparams
      }else{
        if(env2.model==1){
          inparams [c(6,8,14,15,18,31,32,36:39)] <- "NA"
          inputs <- inparams
        }else{
          inparams [c(6,8,14,15,18,31,32,34,35,38,39)] <- "NA"
          inputs <- inparams          
        }
      }
      outputs <- c(excess, critical.res[[3]], critical.res[[4]], "NA", critical.res[[5]], 
                   critical.res[[6]], critical.res[[7]], critical.res[[8]])
   }else{
      mod <- "quantitative"
      if(env1.model==0){
        inparams [c(4,5,7,9,13,16,17,22:26,29)] <- "NA"
        inputs <- inparams
      }else{
        if(env1.model==1){
          inparams [c(4,5,7,9,13,16,17,20,21,25:28)] <- "NA"
          inputs <- inparams
        }else{
          inparams [c(4,5,7,9,13,16,17,20,21,23,24,27,28)] <- "NA"
          inputs <- inparams          
        }
      }
      if(env2.model==0){
        inparams [c(4,5,7,9,13,16,17,33:37,40)] <- "NA"
        inputs <- inparams
      }else{
        if(env2.model==1){
          inparams [c(4,5,7,9,13,16,17,31,32,36:39)] <- "NA"
          inputs <- inparams
        }else{
          inparams [c(4,5,7,9,13,16,17,31,32,34,35,38,39)] <- "NA"
          inputs <- inparams          
        }
      }
      outputs <- c("NA", "NA", "NA", critical.res[[3]], critical.res[[4]], critical.res[[5]], 
                   critical.res[[6]], critical.res[[7]])
   }

   jth.row <- as.character(c(inputs,outputs))
   write(t(jth.row),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
}
}
