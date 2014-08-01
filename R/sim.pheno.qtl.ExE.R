#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param environment1 Exposure data for the first determinant
#' @param env1.efkt Effects of the first determinant
#' @param environment2 Exposure data for the second determinant
#' @param env2.efkt Effects of the second determinant
#' @param interaction data
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.qtl.ExE <-
function(numsubjects=NULL,pheno.mean=NULL,pheno.sd=NULL,environment1=NULL,env1.efkt=NULL,
         environment2=NULL,env2.efkt=NULL,interaction=NULL,int.efkt=NULL)
{  

  # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   num.obs <- numsubjects
   alpha <- pheno.mean 

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (env1.efkt*environment1) + (env2.efkt*environment2) + (int.efkt*interaction)

   # GENERATE THE TRUE PHENOTYPE DATA
   phenotype <- rnorm(num.obs,lp,pheno.sd)
   
   return(phenotype)
}

