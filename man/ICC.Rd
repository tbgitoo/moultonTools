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
\item{method}{Method for calculation of the intraclass correlation coefficient. This can be \code{"unbiased"} (equivalently, \code{"Fisher"}), or \code{"unbiased_unequal_cluster_size"} for
a variant of the Fisher formula that will yield ICC=1 also for unequal cluster sizes upon full clustering, or finally \code{"ANOVA"} to use ANOVA for estimation of the intraclass-correlation.} 
 \item{var_n}{Calculation of the population variance with n or n-1 degrees of freedom. By default,\code{var_n="n"} and n is used. To use n-1, pass \code{var_n="n-1"}.}
 
}


\value{
The intraclass correlation in the \code{outcome} assuming the group structure provided by \code{group}. The intraclass correlation indicates how correlated the values within the groups are relative to the overall variance. The actual values can also be negative (for methods "unbiased", "unbiased_unequal_cluster_size") since negative numbers are needed to obtain an unbiased expectation value of 0. For complete clustering (identical outcome values per group), the "unbiased_unequal_cluster_size" method will produce an ICC value of 1. 
}
\details{In the literature, various calculation methods for the intra-class correlation coefficient. Some are implemented here, and can be chosen by setting the value of the \code{method} argument.\cr\cr

The default value of the \code{method} argument is \code{"unbiased"}, alternatively the same method can be accessed by providing \code{method="Fisher"}. This evaluation indeed corresponds to application of the original Fisher formula. Historically, the initial formula was Fisher's sibling formula (Fisher R.A., Statistical Methods for Research Workers (Twelth Ed.), Edinburgh: Oliver and Boyd), which was then expanded to general cluster size by various authors. The exact formula used here is the formula given as the "px" definition below eq. 8.2.5 in "Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008. This formula is reportedly unbiaised (e.g. Angrist and Pischke) and so it should have an expectation value of 0 for uncorrelated random variables where \code{outcome} has no particular link with the clustering variable \code{group}. One can also verify that it returns 1 for complete clustering (as defined by identical values of \code{outcome} within each level of the \code{group} argument. However, as one also can verify, if the group sizes are not equal, this second property does not hold anymore. \cr\cr

We provide also an unbiased estimator which is equal to exactly equal to 1 for the totally correlated case. It can be invoked by providing \code{method="unbiased_unequal_cluster_size"} to \link{ICC}. This estimator is derived from the ANOVA estimator of intraclass correlation, which relates the intergroup variance to the total variance, but it is corrected for biais to give an expectation value of 0 for independently, equally and normally distributed outcome values. This estimator is \code{(n-3)/(n-N-2)*sum(ng*(zg-zbar))/sum(zi-zbar)-(N-1)/(n-N-2)}. The symboles are: \code{n} is the total number of observations (i.e. \code{length(outcome)=length(group)}), \code{N} the number of groups (i.e. \code{length(unique(group))}), \code{zbar} the overall mean (given by \code{zbar=mean(outcome)}), \code{ng} the individual group sizes (obtainable for example by \code{aggregate(outcome ~ group, FUN=length)}) and \code{zg} are the group means (obtainable for example by \code{aggregate(outcome ~ group, FUN=length)}).

In addition, we provide an ANOVA-type estimator for completeness. For this, \link{ICC} needs to be invoked with \code{method="ANOVA"}. In this case, \link{anova} on a standard linear model \link{lm} is used for variance decomposition, and the fraction of the variance represented by the inter-group variance as compared to the total variance is returned as the intraclass correlation coefficient.\cr\cr

For the methods \code{"unbiased"} and \code{""unbiased_unequal_cluster_size"}, an intergroup correlation term is explicitly compared to the sample variance (see eq. 8.2.5 in "Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008). This raises the question of estimation of the variance, either by \code{sum(squares)/n} or \code{sum(squares)/(n-1}. The estimation of the variances with a denominator of \code{n-1} is more precise on real-world samples, but it causes the estimated intra-class correlation coefficient not to reach is theoretical value of 1 for fully correlated samples. For this reason, by default, \code{n} is used, but we provide the option of forcing the use of a \code{n-1} degrees of liberty for variance estimation by providing \code{var_n="n-1"}. The value of the argument \code{var_n} has no effect when using \code{method="ANOVA"} since then, all relevant sums of squares are calculated using \link{anova} on a linear model as provided by \link{lm}.\cr\cr 

For all methods, the calculation of the intraclass correlation involves sums of squares. No continuity correction is applied, so caution should be used if some of the group sizes are less than 5.

}

\author{
Thomas and Marina Braschler
}
\examples{
# Fully correlated with equal group sizes. For this particular case, one expects the three methods to all give an intraclass correlation coefficient of exactly 1 
outcome=c(0,0,0,1,1,1,-1,-1,-1)
group=c(1,1,1,2,2,2,3,3,3)
ICC(outcome,group) # Fisher formula by default
ICC(outcome,group,method="unbiased_unequal_cluster_size")
ICC(outcome,group,method="ANOVA")
# Fully correlated with unequal group size. In this case, the Fisher formula (by default or by method="unbiased") typically returns results different from 1, whereas the other two methods
# will return precisely 1
outcome=c(0,0,0,1,1,-1,-1,-1,-1,-1,-1)
group=c(1,1,1,2,2,3,3,3,3,3,3)
ICC(outcome,group) # Fisher formula by default
ICC(outcome,group,method="unbiased_unequal_cluster_size")
ICC(outcome,group,method="ANOVA")
# Expectation value (biais) for the uncorrelated case
ICC_values=data.frame("Fisher"=vector(mode="numeric",length=0),"unbiased_unequal_cluster_size"=vector(mode="numeric",length=0),"ANOVA"=vector(mode="numeric",length=0))
for(ind in 1:300)
{
outcome=rnorm(15)
group=c(rep(1,5),rep(2,10))
current=data.frame("Fisher"=vector(mode="numeric",length=1),"unbiased_unequal_cluster_size"=vector(mode="numeric",length=1),"ANOVA"=vector(mode="numeric",length=1))
for(theMethod in colnames(current))
{
	current[,theMethod]=ICC(outcome,group,method=theMethod)
}
ICC_values=rbind(ICC_values,current)
}
lapply(ICC_values,mean)

#If the outcome is a factor, this is assumed unordered and the mean ICC for presence/absence of each level is returned.
#In the complete clustered case, factors levels are distinct and identically associated with the group, the ICC value should be 1
outcome=as.factor(c(rep("A",3),rep("B",2),rep("C",6)))
group=c(1,1,1,2,2,3,3,3,3,3,3)
ICC(outcome,group) # Fisher formula by default
ICC(outcome,group,method="unbiased_unequal_cluster_size")
ICC(outcome,group,method="ANOVA")

\dontrun{
# Further detailed calculations are given in this pdf file within the /inst folder of this package. Here's the path:
system.file("calculations","ICC_calculations.pdf",package="moultonTools")

# To directly view it, this might work on some systems. Otherwise you'll have to navigate to the path given above if interested.
n <- grep("pdf", names(options()))
system_pdf_viewer=options()[[n]]
system_pdf_command = paste(system_pdf_viewer, system.file("calculations","ICC_calculations.pdf",package="moultonTools"),"2> /dev/null &",sep=" ")
system(system_pdf_command)


# Direct comparison to the R-Companion to Angrist and Pischke

library(devtools)
install_github(repo = "MatthieuStigler/RCompAngrist", subdir = "RCompAngrist")
library(RCompAngrist)
outcome=c(-1,-1.5,-2,-1.5,2,2.5,3,2.5,1.5)
group=c(1,1,1,1,2,2,2,2,2)
# This is the ICC function for the uniased method in the Angrist R Companion (it is hidden and has to be accessed using :::)
RcompAngrist:::IC_moult(outcome,group)
# For reference, this is the same formula here, when n-1 is imposed for the evaluation of the variance
ICC(outcome,group,var_n="n-1")



}



}

\seealso{
  The ANOVA approach to intraclass correlation offers substantial flexibility and power to analyse many different situations and models. A simple ANOVA method has been included here for completeness, but it is generally advisable to use some of these more flexible frameworks. They are particularly well established in the field of psychology, as one can see from the available R packages. Some examples are:\cr
  The package "psy" provides an ICC evaluation (function icc in the psy package) for data with identified raters and subjects (here, by equivalence, raters would be known but not the exact subject matching). \cr
  The package "irr" provides ICC evaluation (function icc in the irr package) based on ANOVA, covering many different cases and also hypothesis testing.\cr
  The package "psych" similarly provides extensive ICC evaluation (function ICC in the psych package) based on ANOVA, covering many different cases, definitions and also hypothesis testing.\cr
  The package "rptR" focuses on repeatibility aspects and adds permutation tests and confidence intervals.\cr
  \cr\cr
  More specifically regarding econometrics, there is a R companion to the book "Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008. It is available on as R package on Github (see https://github.com/MatthieuStigler/RCompAngrist), although it seems not to have been maintained further.
 
}

\references{

"Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008 \cr\cr

"Statistical Methods for Research Workers", Fisher R.A.,(Twelth Ed.), Edinburgh: Oliver and Boyd, 1954


}

\keyword{ misc }

