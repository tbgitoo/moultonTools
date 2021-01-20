\name{anova_moulton}
\alias{anova_moulton}

\title{
anova_moulton
}
\description{
Anova for analysis including clustered variables. Corrects the F values by taking into account intraclass correlation for known grouping by the Moulton formula (Mostly Harmless Econometrics: An Empiricist's Compagnion, Joshua D. Angrist and Joern-Steffen Pischke), and corrects the degrees of freedom by the between-within method (Singer JD. Using SAS PROC MIXED to fit multilevel models, hierarchical models, and individual growth models. Journal of Educational and Behavioral Statistics. 1998;24:323-355).
}
\usage{
anova_moulton(linmod,data_object,grouping_col,use_predictors=TRUE,method="unbiased")
}

\arguments{
  \item{linmod}{
	Standard linear model based on \link{lm}. Do not include the clustering variable in the model
}

  \item{data_object}{
  The underlying data that is also used in \code{linmod}
}
 \item{grouping_col}{The name of the column in \code{data_object} that indicates the clustering. Thus, \code{data_object[,"grouping_col"]} should be a factor indicating which observation belongs to which cluster.}
\item{use_predictors}{Switch indicating whether the linear predictors in \code{linmod} should also be used in the calculation of the Moulton factor. In general, they should,  \code{use_predictors} should be use for debugging purposes only}
\item{method}{Method argument for \link{moulton_factor} estimation, ultimately passed on to \link{ICC}}
}


\value{
An object of class ANOVA, including F and P values adjusted for the clustering by means of the Moulton method and between-within adjustment of degrees of freedom. The moulton factors themselves are returned as \code{attr(x,"moulton")}, whereas the degrees of freedom are returned as attributes \code{attr(x,"df2")}, assuming \code{x} to be the output variable.
    
    
}

\author{
Marina + Thomas Braschler
}

\examples{

# Typically example of clustering. If one looks only at the linear regression between the outcome and regressor columns (data_object below), there is a clear trend and also a relatively
# large number of datapoints. As a consequence, simple linear regression gives a significant result
data_object = data.frame(outcome=c(1,2,5,6,12,5,2,3,20,19,10,15,13,20,14),regressor=c(1,1,1,3,3,3,7,7,7,7,10,10,10,10,10),clusters=c("A","A","B","A","B","B","B","B","C","C","C","C","C","C","C"))
linmod=lm(outcome~regressor,data_object)
summary(linmod)
# However, knowlege of the cluster column (last column in the data_object variable, labelled "clusters") indicates that this might also be because the regressor values 
# tend to be associated with the clusters "A"-"C"
# Intraclass correlation analysis indeed shows that there is a tendency for clustering of the outcome data along the cluster variable
ICC(data_object$outcome,data_object$clusters)
# and also clustering of the regressor values within the clusters
ICC(data_object$regressor,data_object$clusters)

# This however means that outcome values are partial replicates per cluster, which in addition line up unfortunately with the regressors because of the very imbalanced experimental design. 
# And indeed, Moulton analysis shows that the result
# should probably not be considered significant
anova_moulton(linmod,data_object,"clusters")
# One gets to a similar conclusion by aggregating both the outcome and the regressors within the clusters: The trend seems to be there, but the effective n is too low
aggregated_data_object=aggregate(data_object[,c("outcome","regressor")],by=list(clusters=data_object$clusters),FUN=mean)
summary(lm(outcome~regressor,aggregated_data_object))

}


\references{

"Mostly harmless econometrics: An Empiricist's Companion", Angrist J.D., Pischke, J.S., 2008 \cr\cr

"Using SAS PROC MIXED to fit multilevel models, hierarchical models, and individual growth models." Singer JD.  Journal of Educational and Behavioral Statistics. 1998;24:323-355

}

\keyword{ misc }

