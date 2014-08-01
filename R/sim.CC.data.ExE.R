#' 
#' @title Simulates case and controls
#' @description Generates affected and non-affected subjects until the set sample 
#' size is achieved.
#' @param n Number of observations to generate per iteration
#' @param ncases Number of cases to simulate
#' @param ncontrols Number of controls to simulate
#' @param max.sample.size Maximum number of observations allowed
#' @param pheno.prev Prevalence of the binary outcome
#' @param e1.model Model of the 1st environmental exposure
#' @param e1.prev Prevelance of the 1st environmental exposure
#' @param e1.mean Mean of the 1st environmental exposure
#' @param e1.sd Standard deviation of the 1st environmental exposure
#' @param e1.low.lim Lower limit of the 1st environmental exposure
#' @param e1.up.lim Upper limit of the 1st environmental exposure
#' @param e1.OR Odds ratios of the 1st environmental exposure
#' @param e2.model Model of the 2nd environmental exposure
#' @param e2.prev Prevelance of the 2nd environmental exposure
#' @param e2.mean Mean of the 2nd environmental exposure
#' @param e2.sd Standard deviation of the 2nd environmental exposure
#' @param e2.low.lim Lower limit of the 2nd environmental exposure
#' @param e2.up.lim Upper limit of the 2nd environmental exposure
#' @param e2.OR Odds ratios of the 2nd environmental exposure
#' @param i.OR Odds ration of the interaction
#' @param b.OR Baseline odds ratio for subject on 95 percent population 
#' centile versus 5 percentile. This parameter reflects the heterogeneity in disease 
#' risk arising from determinates that have not been measured or have not been 
#' included in the model.
#' @param ph.error misclassification rates: 1-sensitivity and 1-specificity
#' @return A matrix
#' @keywords internal
#' @author Gaye A.
#'
sim.CC.data.ExE <-
function(n=NULL, ncases=NULL, ncontrols=NULL, max.sample.size=NULL, pheno.prev=NULL,
e1.model=NULL, e1.prev=NULL, e1.mean=NULL, e1.sd=NULL, e1.low.lim=NULL, e1.up.lim=NULL, 
e1.OR=NULL, e2.model=NULL, e2.prev=NULL, e2.mean=NULL, e2.sd=NULL, e2.low.lim=NULL, 
e2.up.lim=NULL, e2.OR=NULL,i.OR=NULL, b.OR=NULL, ph.error=NULL)
{
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX TO STORE THE GENERATED DATA
   sim.matrix <- matrix(numeric(0), ncol=4)

   # SET LOOP COUNTER
   numloops <- 0

   # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
   # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
   while(complete==0 && complete.absolute==0)
     {

       # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA FOR THE FIRST DETERMINANT
       env1 <- sim.env.data(num.obs=n, env.model=e1.model, env.prev=e1.prev, env.mean=e1.mean, 
                                   env.sd=e1.sd, env.low.lim=e1.low.lim,env.up.lim=e1.up.lim)
       
       # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA FOR THE SECOND DETERMINANT
       env2 <- sim.env.data(num.obs=n, env.model=e2.model, env.prev=e2.prev, env.mean=e2.mean, 
                           env.sd=e2.sd, env.low.lim=e2.low.lim, env.up.lim=e2.up.lim)
       
       # GENERATE THE TRUE INTERACTION DATA
       int <- env1 * env2       

       # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINE RISK: 
       # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
       # VARIANCE ON SCALE OF LOG-ODDS
       s.effect.data <- sim.subject.data(n, b.OR)
					  
       # GENERATE THE TRUE OUTCOME DATA
       pheno.data <- sim.pheno.bin.ExE(num.obs=n, disease.prev=pheno.prev, environment1=env1, environment2=env2, 
                                       interaction=int, subject.effect.data=s.effect.data, env1.OR=e1.OR, 
                                       env2.OR=e2.OR, int.OR=i.OR)
       true.phenotype <- pheno.data
       
       # GENERATE THE OBSERVED OUTCOME DATA FROM WHICH WE SELECT CASES AND CONTROLS
       obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=0, pheno.error=ph.error)
       pheno <- obs.phenotype
       
       # STORE THE TRUE OUTCOME, GENETIC AND ENVIRONMENT AND ALLELE DATA IN AN OUTPUT MATRIX 
       # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
       sim.matrix.temp <- cbind(pheno,env1,env2,int)

       # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
       sim.matrix <- rbind(sim.matrix, sim.matrix.temp)

       # SELECT OUT CASES
       sim.matrix.cases <- sim.matrix[pheno==1,]

       # SELECT OUT CONTROLS
       sim.matrix.controls <- sim.matrix[pheno==0,]

       # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
       cases.simulated <- dim(sim.matrix.cases)[1]
       controls.simulated <- dim(sim.matrix.controls)[1]

       # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
       if(cases.simulated >= ncases)
       {
         sim.matrix.cases <- sim.matrix.cases[1:ncases,]
         cases.complete <- 1
       }

       # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
       if(controls.simulated>=ncontrols)
       {
         sim.matrix.controls <- sim.matrix.controls[1:ncontrols,]
         controls.complete <- 1
       }

       # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
       complete <- cases.complete*controls.complete		

       # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
       complete.absolute <- (((block+1)*n)>=max.sample.size)
       if(complete.absolute==1) {sample.size.excess <- 1}else{sample.size.excess <- 0}

        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
       
   }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)
   # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "environment1", "environment2", "interaction")
   mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess)
}

