\name{moulton_factor}
\alias{moulton_factor}

\title{
Moulton factor for correction of intraclass correlation in linear regression
}
\description{
The moulton factor corrects for intraclass correlation in linear regression in order not to underestimate the standard deviations of the estimators
}
\usage{
moulton_factor(outcome,group)
}
\arguments{
  \item{outcome}{
 The outcome variable (numerical vector) }
  \item{group}{
The grouping variable, typically a factor indicating to which group the outcomes belong }
 \item{estimator}{Optionally, an estimator (predictive variable in the linear regression) can also be provided along with the outcome and group to get a more realistic picture of the Moulton factor 
}
}

\value{
A number generally larger than 1, indicating the correction factor for the standard deviations of the estimated slopes arising from the grouping structure
}
\details{The package uses  \link{ICC} internally }
\author{
Thomas and Marina Braschler
}
\examples{
data(acuity_data_extended) 
moulton_factor(acuity_data_extended$logMAR[acuity_data_extended$month<=48],acuity_data_extended$patient[acuity_data_extended$month<=48])

outcome=acuity_data_extended$logMARNormalized[acuity_data_extended$month==12 & acuity_data_extended$eye != "ND"]

group=acuity_data_extended$eye[acuity_data_extended$month==12 & acuity_data_extended$eye != "ND"]
moulton_factor(outcome,group)
}

\references{Mostly Harmless Econometrics: An Empiricist's Compagnion, Joshua D. Angrist and Joern-Steffen Pischke}

\keyword{ misc }

