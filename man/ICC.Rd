\name{ICC}
\alias{ICC}

\title{
Intraclass-correlation ICC
}
\description{
Estimates the intraclass correlation coefficient for a given outome (readout) and grouping (clustering) variable
}
\usage{
ICC(outcome,group,method="unbiased")
}
\arguments{
  \item{outcome}{
 The outcome variable (numerical vector) }
  \item{group}{
The grouping variable, typically a factor indicating to which group the outcomes belong }
 
}


\value{
A number with expections values between 0 and 1, indicating how correlated the values within the groups are. The actual values can also be negative (for methods "unbiased", "unbiased_unequal_cluster_size") since negative numbers are needed to obtain an unbiased expectation value of 0. For complete clustering (identical outcome values per group), the "unbiased_unequal_cluster_size" method will produce an ICC value of 1. 
}
\details{In the literature, various calculation methods for the intra-class correlation coefficient. Some are implemented here, and can be chosen by setting the value of the \code{method} argument.\cr\cr

The default value of the \code{method} argument is "unbiased". This corresponds to application of the original Fisher formula: Historically, the initial formula was Fisher's sibling formula (Fisher R.A., Statistical Methods for Research Workers (Twelth Ed.), Edinburgh: Oliver and Boyd), which was then expanded to general cluster size by various authors. The exact formula used here is the formula given as the "px" definition below eq. 8.2.5 in "Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008. This formula is reportedly unbiaised (e.g. Angrist and Pischke) and so it should have an expectation value of 0 for uncorrelated random variables where \code{outcome} has no particular link with the clustering variable \code{group}. One can also verify that it returns 1 for complete clustering (as defined by identical values of \code{outcome} within each level of the \code{group} argument. However, as one also can verify, if the group sizes are not equal, this second property does not hold anymore. \cr\cr

 in The package uses the \link{multilevel} regression function }
\author{
Thomas and Marina Braschler
}
\examples{
data(acuity_data_extended) 
ICC(acuity_data_extended$logMAR[acuity_data_extended$month<=48],acuity_data_extended$patient[acuity_data_extended$month<=48])

outcome=acuity_data_extended$logMARNormalized[acuity_data_extended$month==12 & acuity_data_extended$eye != "ND"]

group=acuity_data_extended$eye[acuity_data_extended$month==12 & acuity_data_extended$eye != "ND"]
ICC(outcome,group)
}

\keyword{ misc }

