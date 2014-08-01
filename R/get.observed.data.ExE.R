#'
#' @title Generates exposure data with some error
#' @description Uses functions make.obs.geno and make.obs.env to generate effect data with a set level of error
#' @param data Input table of simulated data considered as true data
#' @param e1.error Misclassification rates in the 1st environmental exposure assessment: 1-sensitivity and 1-specificity
#' @param e1.model Model of the 1st exposure: binary=0, quantitative-normal=1 or quantitative-uniform=2
#' @param e1.prev Prevalence of the 1st environmental exposure
#' @param e1.sd Standard Deviation of the 1st environmental exposure
#' @param e1.reliability Reliability of the assessment of 1st environmental exposure
#' @param e2.error Misclassification rates in the 2nd environmental exposure assessment: 1-sensitivity and 1-specificity
#' @param e2.model Model of the 2nd exposure: binary=0, quantitative-normal=1 or quantitative-uniform=2
#' @param e2.prev Prevalence of the 2nd environmental exposure
#' @param e2.sd Standard Deviation of the 2nd environmental exposure
#' @param e2.reliability Reliability of the assessment of 2nd environmental exposure
#' @return A matrix
#' @keywords internal
#' @author Gaye A
#'
get.observed.data.ExE <-
function(data=NULL,e1.error=NULL,e1.model=NULL,e1.prev=NULL,e1.sd=NULL, e1.reliability=NULL,
         e2.error=NULL,e2.model=NULL,e2.prev=NULL,e2.sd=NULL, e2.reliability=NULL)
{
		sim.df <- data      

    # GET THE OBSERVED ENVIRONMENTAL EXPOSURE DATA FOR THE FIRST DETERMINANT
    true.environment1 <- sim.df$environment1
    obs.environment1 <- get.obs.env(env.data=true.environment1,env.model=e1.model,env.prev=e1.prev,
                                   env.sd=e1.sd,env.error=e1.error,env.reliability=e1.reliability)
    
		# GET THE OBSERVED ENVIRONMENTAL EXPOSURE DATA FOR THE SECOND DETERMINANT
		true.environment2 <- sim.df$environment2
		obs.environment2 <- get.obs.env(env.data=true.environment2,env.model=e2.model,env.prev=e2.prev,
		                                env.sd=e2.sd,env.error=e2.error,env.reliability=e2.reliability)
    
    # GET THE OBSERVED INTERACTION DATA
    obs.interaction <- obs.environment1 * obs.environment2

    # REPLACE THE TRUE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
    # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
    sim.df$environment1 <- obs.environment1
		sim.df$environment2 <- obs.environment2
		sim.df$interaction <- obs.interaction
    
		# RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    colnames(sim.df) <- c("id", "phenotype", "environment1", "environment2", "interaction")
    return(sim.df)
}

