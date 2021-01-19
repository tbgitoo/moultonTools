\name{t.test.moulton}
\alias{t.test.moulton}

\title{
t.test.moulton
}
\description{
T-test with Moulton correction for clustered variables. Works like \link{t.test}, but corrects the standard deviation and in some cases the degrees of freedom applicable in test evaluation
}
\usage{
t.test.moulton(x, y = NULL,cluster_x=NULL,cluster_y=NULL, alternative = c("two.sided", "less", "greater"),
mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
...)
}

\arguments{
  \item{x}{
	a (non-empty) numeric vector of data values
}
\item{y}{
	an optional (non-empty) numeric vector of data values
}

 \item{cluster_x}{
	Cluster identities for the values in  x
}
 \item{cluster_y}{
	Cluster identities for the values in  y
}
\item{alternative}{
	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test).}
\item{paired}{a logical indicating whether you want a paired t-test}
\item{var.equal}{a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.}
\item{conf.level}{confidence level of the interval.}

\item{...}{If NULL is provided explicitly for both \code{cluster_x} and \code{cluster_y}, the \link{t.test} function is called and the arguments \code{...} are passed on to \link{t.test} }

}
\details{ There are a number of different cases for the t-test, the following corrections are applied: \cr\cr
(1) If \code{y=NULL}, the variable \code{x} is tested against the theoretical mean \code{mu}. In that case, the Moulton correction factor is obtained from the values of \code{x} and the clustering of \code{x} given by \code{cluster_x}. The degree of freedom is estimated as the number of clusters minus 1 (i.e. \code{length(unique(cluster_x))-1}.\cr
(2) If both \code{x} and \code{y} are provided, the further approach depends on whether a paired or unpaired test is desired. This depends on the value of the \code{paired} argument. If a paired test is desired, \code{x} is replaced by \code{x-y}, and only the clustering described in \code{cluster_x} is used, the approach being otherwise identical to the case (1) given above. \cr
(3) If both \code{x} and \code{y} are provided, and an unpaired test is desired \code{paired=FALSE}, then the clustering of \code{x} and \code{y} is assumed to be described separately by \code{cluster_x} and \code{cluster_y}. In fact, in this case, \code{cluster_x} and \code{cluster_y} must no overlap, because the clustering is supposed to take place within the \code{x} and \code{y}, not in a crossed manner between them. In the case, the Mouton correction is applied to the joint variance; the degree of freedom is corrected according to the number of clusters rather than the number of observations; the details depend on whether \code{var.equal} is provided \code{TRUE} or \code{FALSE} (see \link{t.test} for details about this)   }

\value{
An object of class \code{"htest"}, see \link{t.test}. In addition, 
    
    
}

\author{
Marina + Thomas Braschler
}

\examples{


# Example 1: Strong clustering vs. no clustering


# x and y random variables

n=10


cluster_sd = 1

cluster_values = rnorm(6,cluster_sd)

cluster_x=c(rep(1,n),rep(2,n),rep(3,n))
cluster_y=c(rep(4,n),rep(5,n),rep(6,n))
x=cluster_values[cluster_x]+rnorm(3*n,sd=0.15)
y=cluster_values[cluster_y]+rnorm(3*n,sd=0.15)

y=c(cluster_values[4]+rnorm(n,sd=0.15),cluster_values[5]+rnorm(n,sd=0.15),cluster_values[6]+rnorm(n,sd=0.15))


# Simple t-testing gives quite often a highly significant result

t.test(x,y)

# If we know that x and y come each from only four distinct clusters, then the p-values should actually be closer to what we observe at the level of the clusters:

# While really, we should have sampled on the level of the clusters
t.test(aggregate(x ~ cluster_x,FUN=mean)$x,aggregate(y ~ cluster_y,FUN=mean)$y)

# The Moulton correction roughly does it

t.test.moulton(x,y,cluster_x,cluster_y)

# If on the other hand, the groups are distributed over the two test variables, this is a balanced design and chances to get a highly significant result where there is nothing are not that high anymore

cluster_x=c(rep(1,n),rep(2,n),rep(3,n))
cluster_y=cluster_x
x=cluster_values[cluster_x]+rnorm(3*n,sd=0.15)
y=cluster_values[cluster_y]+rnorm(3*n,sd=0.15)



t.test(x,y)

t.test.moulton(x,y,cluster_x,cluster_y)







}

\keyword{ misc }

