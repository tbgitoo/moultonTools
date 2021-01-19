\name{ICC}
\alias{ICC}

\title{
Intraclass-correlation ICC
}
\description{
Estimates the intraclass correlation using no other predictor than the group identity in linear regression
}
\usage{
ICC(outcome,group)
}
\arguments{
  \item{outcome}{
 The outcome variable (numerical vector) }
  \item{group}{
The grouping variable, typically a factor indicating to which group the outcomes belong }
 
}


\value{
A number between 0 and 1, indicating how correlated the values within the groups are
}
\details{The package uses the \link{multilevel} regression function }
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

