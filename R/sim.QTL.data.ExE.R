#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param e1.model Model of the 1st environmental exposure
#' @param e1.efkt Effects of the 1st environmental exposure
#' @param e1.prev Prevalence of the 1st environmental exposure
#' @param e1.mean Mean of the 1st environmental exposure
#' @param e1.sd Standard deviation of the 1st environmental exposure
#' @param e1.low.lim lower limit of the 1st environmental exposure
#' @param e1.up.lim upper limit of the 1st environmental exposure
#' @param e2.model Model of the 2nd environmental exposure
#' @param e2.efkt Effects of the 2nd environmental exposure
#' @param e2.prev Prevalence of the 2nd environmental exposure
#' @param e2.mean Mean of the 2nd environmental exposure
#' @param e2.sd Standard deviation of the 2nd environmental exposure
#' @param e2.low.lim lower limit of the 2nd environmental exposure
#' @param e2.up.lim upper limit of the 2nd environmental exposure
#' @param i.efkt Interaction effect
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.

sim.QTL.data.ExE <-
function(n=NULL,ph.mean=NULL,ph.sd=NULL,e1.model=NULL,e1.efkt=NULL,e1.prev=NULL,e1.mean=NULL,e1.sd=NULL,
         e1.low.lim=NULL,e1.up.lim=NULL,e2.model=NULL,e2.efkt=NULL,e2.prev=NULL,e2.mean=NULL,e2.sd=NULL,
         e2.low.lim=NULL,e2.up.lim=NULL,i.efkt=NULL,
         pheno.rel=NULL)
{
	   
   # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA FOR THE FIRST DETERMINANT			
   env1 <- sim.env.data(num.obs=n,env.model=e1.model,env.prev=e1.prev,env.mean=e1.mean,
                       env.sd=e1.sd,env.low.lim=e1.low.lim,env.up.lim=e1.up.lim)
   
   # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA FOR THE SECOND DETERMINANT    	
   env2 <- sim.env.data(num.obs=n,env.model=e2.model,env.prev=e2.prev,env.mean=e2.mean,
                        env.sd=e2.sd,env.low.lim=e2.low.lim,env.up.lim=e2.up.lim)  
   
   # GENERATE THE TRUE INTERACTION DATA           
   int <- env1 * env2
           
   # GENERATE THE TRUE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.ExE(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,
                                   environment1=env1,env1.efkt=e1.efkt,
                                   environment2=env2,env2.efkt=e2.efkt,
                                   interaction=int,int.efkt=i.efkt)
   true.phenotype <- pheno.data
   
   # GENERATE THE OBSERVED OUTCOME DATA 
   obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                  pheno.sd=ph.sd, pheno.reliability=pheno.rel)
   pheno <- obs.phenotype
   

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(pheno,env1,env2,int)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype", "environment1","environment2","interaction")
   mm <- data.frame(sim.matrix)
}

