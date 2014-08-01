#'
#' @title Generates phenotype status
#' @description Generates affected and non-affected subjects
#' @param num.obs Number of observations to generate per iteration
#' @param disease.prev Prevalence of the binary outcome
#' @param environment1 Exposure data for the 1st determinant
#' @param environment2 Exposure data for the 2nd determinant
#' @param interaction data
#' @param subject.effect.data Subject effect data, reflects the heterogenity 
#' in baseline disease risk
#' @param env1.OR Odds ratios of the 1st determinant
#' @param env2.OR Odds ratios of the 2nd determinant
#' @param int.OR Odds ration of the interaction
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.bin.ExE <-
function(num.obs=NULL, disease.prev=NULL, environment1=NULL, environment2=NULL, interaction=NULL, 
         subject.effect.data=NULL, env1.OR=NULL, env2.OR=NULL, int.OR=NULL)
{ 
   # GET THE ALPHA AND BETA VALUES
   alpha <- log(disease.prev/(1-disease.prev))
   env1.beta <-  log(env1.OR)   
   env2.beta <-	log(env2.OR)
   int.beta <- log(int.OR)

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (env1.beta*environment1) + (env2.beta*environment2) + (int.beta*interaction) + subject.effect.data
   # GET THE 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
   mu <- exp(lp)/(1 + exp(lp))
   
   # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
   phenotype <- rbinom(num.obs,1,mu)
   
   return(phenotype)
}

