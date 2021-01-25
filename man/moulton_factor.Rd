\name{moulton_factor}
\alias{moulton_factor}

\title{
Moulton factor for correction of intraclass correlation in linear regression
}
\description{
The moulton factor corrects for intraclass correlation in linear regression in order not to underestimate the standard deviations of the estimators
}
\usage{
moulton_factor(outcome,group,estimator=NULL,...)
}
\arguments{
  \item{outcome}{
 The outcome variable (numerical vector) }
  \item{group}{
The grouping variable, typically a factor indicating to which group the outcomes belong }
 \item{estimator}{Optionally, an estimator (predictive variable in the linear regression) can also be provided along with the outcome and group to get a more realistic picture of the Moulton factor }
 \item{...}{Additional parameters to be passed to \link{ICC}, called internally. Typically, this is a \code{method} argument as documented for the intraclass correlation coefficient \link{ICC}, possible values are at present \code{method="unbiased"} or equivalently \code{method="Fisher"}, as well as \code{method="unbiased_unequal_cluster_size"} and \code{"ANOVA"}. It is also possible to provide \code{method="Angrist_mixed"}. In this case, the approach suggested in "Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008 is, where for the intraclass correlation of the observed values (i.e. \code{outcome}), the ANOVA approach is used, whereas for the regressors (i.e. \code{estimators}), Fisher's approach is used.
}
}

\value{
A number equal to or larger than 1, indicating the correction factor for the standard deviations of the estimated slopes arising from the grouping structure. Values smaller than 1 (arising from negative intraclass correlation estimates) are set equal to 1 with the aim to maintain test statistics neutral or conservative in all cases. The intraclass correlation coefficient of the \code{outcome} variable under the cluster structure given by \code{group} is returned as attribute "ICC" to the primary return value, while the intraclass correlation coefficient for the \code{estimator} variable, also under the grouping structure given by \code{group} is returned as attribute "px".
}
\details{The package uses  \link{ICC} internally. Of note, balancing designs regarding regression variables are important, and if efforts were made to obtain balanced designs, the regressor values should be supplied as \code{estimator} arguments. If no \code{estimator} argument is supplied, the worst-case scenario is assumed where regressors are entirely specific to the clusters.\cr\cr
Also, Mouton factors can be used to correct linear statistics such as t-statistics, or squared statistics, such as chisquared values. The value provided here is for correction of the linear statistics; if application to Chisquare variables is desired, the square needs to be taken. }
\author{
Thomas and Marina Braschler
}

\seealso{
The Moulton function draws on the functionality of intraclass-correlation, see \link{ICC}.
}

\examples{
# Worst case scenario with a priori unknown or badly balanced regressors
outcome=c(-1,-1.5,-2,-1.5,2,2.5,3,2.5)
group=c(1,1,1,1,2,2,2,2)
moulton_factor(outcome,group)
# If a balanced design is applied for a regression (or, with only two levels, equivalently a t-test), there is no clustering issue and a Moulton factor close to 1 should result:
outcome=c(-1,-1.5,-2,-1.5,2,2.5,3,2.5)
group=c(1,1,1,1,2,2,2,2)
estimator=c(0,0,1,1,0,0,1,1)
moulton_factor(outcome,group,estimator=estimator)
# Supply a method argument to use a specific method for intraclass correlation estimation
outcome=c(-1,-1.5,-2,-1.5,2,2.5,3,2.5,1.5)
group=c(1,1,1,1,2,2,2,2,2)
estimator=c(0,0,1,1,0,0,1,1,1)
moulton_factor(outcome,group,estimator=estimator,method="ANOVA")
# Specifically for ANOVA, it can be useful if the estimator (aka, regressors in this case) are factors. As in:
outcome=c(-1,-1.5,-2,-1.5,2,2.5,3,2.5,1.5)
group=c(1,1,1,1,2,2,2,2,2)
estimator=as.factor(c("A","A","B","B","C","C","A","A","A"))
moulton_factor(outcome,group,estimator=estimator)




}

\references{"Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008}

\keyword{ misc }

