\name{prop.trend.test.moulton}
\alias{prop.trend.test.moulton}

\title{
prop.trend.test.moulton
}
\description{
Proportion trend test, with Moulton correction for clustering. 
}
\usage{
prop.trend.test.moulton(binaryResponse, score, cluster=seq_along(binaryResponse) )
}

\arguments{
  \item{binaryResponse}{
	A vector of \code{TRUE} and \code{FALSE} values ("logical")
}
  \item{score}{
	A numerical vector of the same length as \code{binaryResponse}. This variable is used to evaluate the proportion of succesful (\code{TRUE}) outcomes for each value of \code{score}; these proportions are the basis for \link{prop.trend.test} invoked by this function. 
}

\item{cluster}{
	A vector indicating the cluster identities, of the same length as \code{binaryResponse}
}
\item{...}{Additional arguments to be passed to the underlying \link{moulton_factor} function, particularly a \code{method} argument to use different methods for the estimation of intraclass correlation. See also \link{ICC}}.

}
\details{ The function evaluates the proportions from the binary response variable (\code{binaryResponse}) categorized using the \code{score} variable. It then invokes \link{prop.trend.test} with these proportions. Then, the Moulton correction is applied, by means of a Moulton factor calculated by \link{moulton_factor} based on the \code{cluster} column }

\value{
An object of class \code{"htest"}, see \link{prop.trend.test}. In addition, the Moulton factor used is indicated in a "moulton_factor" field.
    
    
}

\author{
Marina + Thomas Braschler
}

\examples{


# Example 1: Construction of clustered variables and testing


# The primary question about Moulton correction is false positive detection

\dontrun{

one_round_simulation<-function(p0=0.3,sd_cluster=1,n_p=3,n_cluster_per_p=2,n_per_cluster=10)
{


x0 = log(p0/(1+p0))  # To scale the interval 0,1 to -inf to +inf


xtot = vector(mode="numeric",length=0)
score =  vector(mode="numeric",length=0)
cluster_id =  vector(mode="numeric",length=0)


for(ind_score in 1:n_p)
{


score = c(score, rep(ind_score,n_per_cluster*n_cluster_per_p))

	for(ind_cluster in 1:n_cluster_per_p)
	{
		cluster_id = c(cluster_id, rep((ind_score-1)*n_cluster_per_p+ind_cluster,n_per_cluster))
		xtot = c(xtot,rep(x0+rnorm(1,sd=sd_cluster),n_per_cluster))
	}

}

ptot=exp(xtot)/(1+exp(xtot))


response = vector(mode="numeric",length=0)


for(ind in 1:length(ptot))
{
	response=c(response,rbinom(1,1,ptot[ind]))
}


x_prop_trend_test = aggregate(response ~ score, FUN=sum)$response
n_prop_trend_test = aggregate(response ~ score, FUN=length)$response

p_vector=vector(mode="numeric",length=2)

names(p_vector)=c("no.correction","moulton")

p_vector["no.correction"]=prop.trend.test(x_prop_trend_test,n_prop_trend_test)$p.value

p_vector["moulton"]=prop.trend.test.moulton(response, score, cluster_id )$p.value

return(p_vector)


}

p_matrix=matrix(nrow=500,ncol=2)

colnames(p_matrix)=c("no.correction","moulton")


for(ind in 1:(dim(p_matrix)[1]))
{

	p_matrix[ind,]=one_round_simulation(p0=0.3,sd_cluster=0.5,n_p=
	3,n_cluster_per_p=2,n_per_cluster=40)


}



cat(paste("Fraction of false positives at alpha = 0.05 WITHOUT Moulton correction: ",sum(p_matrix[,1]<0.05)/(dim(p_matrix)[1]),sep=""))

cat(paste("Fraction of false positives at alpha = 0.05 WITH Moulton correction: ",sum(p_matrix[,2]<0.05)/(dim(p_matrix)[1]),sep=""))

}

}

\keyword{ misc }

