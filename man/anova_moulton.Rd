\name{anova_moulton}
\alias{anova_moulton}

\title{
anova_moulton
}
\description{
Anova for analysis including clustered variables. Corrects the F values by taking into account intraclass correlation for known grouping by the Moulton formula (Mostly Harmless Econometrics: An Empiricist's Compagnion, Joshua D. Angrist and Joern-Steffen Pischke), and corrects the degrees of freedom by the between-within method (Singer JD. Using SAS PROC MIXED to fit multilevel models, hierarchical models, and individual growth models. Journal of Educational and Behavioral Statistics. 1998;24:323?355).
}
\usage{
anova_moulton(linmod,data_object,grouping_col,use_predictors=TRUE)
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
}
\details{}

\value{
An object of class ANOVA, including F and P values adjusted for the clustering by means of the Moulton method and between-within adjustment of degrees of freedom. The moulton factors themselves are returned as \code{attr(x,"moulton")}, whereas the degrees of freedom are returned as attributes \code{attr(x,"df2")}, assuming \code{x} to be the output variable.
    
    
}

\author{
Marina + Thomas Braschler
}

\examples{

periods_boundaries = c(15,51)
anv = anova_moulton(periods_boundaries=periods_boundaries, combine_HTA_DM=FALSE)
anv

}

\keyword{ misc }

