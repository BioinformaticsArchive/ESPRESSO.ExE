\name{run.espresso.ExE}
\alias{run.espresso.ExE}
\title{Runs a full ESPRESSO analysis}
\usage{
run.espresso.ExE(simulation.params = NULL, pheno.params = NULL,
  env1.params = NULL, env2.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{env1.params}{parameters for the 1st environmental
  determinant}

  \item{env2.params}{parameters for the 2nd environmental
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wish to analyse}
}
\value{
a summary table that contains both the input parameters and
the results of the analysis
}
\description{
This function calls the functions required to run a full
ESPRESSO analysis where the model consists of an outcome
(binary or continuous) determined by two interacting
covariates (a SNP and an environmental exposure)
}
\details{

}
\examples{
{

# load the table that hold the input parameters; each of the table
# hold parameters for 4 scenarios:
# scenario 1: a binary outcome determined by two binary environmental exposures
# and an interaction between the two
# scenario 2: a binary outcome determined by two continuous environmental exposures
# and an interaction between the two
# scenario 3: a quantitative outcome determined by two binary environmental exposures
# and an interaction between the two
# scenario 4: a quantitative outcome determined by two continuous environmental exposures
# and an interaction between the two
data(simulation.params)
data(pheno.params)
data(env1.params)
data(env2.params)

# run the function for the first two scenarios, two binomial models
run.espresso.ExE(simulation.params, pheno.params, env1.params, env2.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.ExE(simulation.params, pheno.params, env1.params, env2.params, scenarios2run=c(3,4))

}
}
\author{
Gaye A.
}

