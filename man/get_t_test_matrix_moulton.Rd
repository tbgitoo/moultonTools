\name{get_t_test_matrix_moulton}
\alias{get_t_test_matrix_moulton}

\title{
Get a T-Test matrix comparing the results under different treatments, taking into account a clustering factor via the Moutlon correction
}
\description{
Compares the data obtained for each value of the controlling factor, and constructs a matrix with the t-test results. The tests are corrected for clustering via the Moulton correction
}
\usage{
get_t_test_matrix_moulton(treatment_factor, data,cluster_factor, ...)
}

\arguments{
  \item{treatment_factor}{
  Factor indicating the treatment conditions (vector)
}
  \item{data}{
 Data obtained under the different conditions (vector of the same length as \code{treatment_factor}
}
\item{cluster_factor}{Factor describing the clustering. This is a vector of the same length as \code{treatment_factor} and \code{data}. Attention: clustering can either concern values inside treatment groups only or clusters across treatment groups as well. The two cases are distinguished by the values passed in \code{cluster_factor}. If they are distinct between two groups (i.e. if the first treatment group has the clusters 1,2,3 and the second the clusters 4,5,6), then the analysis considers clustering separately in the groups; if they are identical between the groups (i.e. clusters 1,2,3 for group 1 and group 2), then the analysis takes the entire clustering structure into account. Mixed scenarios are possible as well (i.e. clusters 1,2,3 in group 1, clusters 1,2 and 4 in group 2).}
\item{...}{
Additional arguments to be passed down to the underlying \code{\link{moulton.t.test}} function
}
}
\details{
The full list of data should be given, so that standard deviations can be calculated. More precisely, there should be several entries where \code{treatment_factor} has identical values. Care should be taken in the description of the \code{cluster_factor}, namely one should avoid to inadvertently imply clusters across groups if this was not the case in the measurement scheme (see the arguments section); this would in particular happen if one was to restart numbering the clusters in each treatment group, while actually the clusters are truly different experimental conditions. 
}
\value{
A matrix containing the p-values of t-tests, comparing the data for each of the possible pairwise combinations of \code{treatment_factor}
}
\author{
Thomas Braschler
}

\seealso{
This is a convenience interface to \link{moulton.t.test}, see \link{moulton.t.test} for further documentation on the actual testing.
}

\examples{


# First, the indicator variables treatment_factor and cluster_factor.
# Three distinct treatments
treatment_factor = c(rep("T1",20),rep("T2",20),rep("T3",20))

# A total of 12 different clusters, 4 per treatment, each with 5 samples
cluster_factor = c(rep("C1",5),rep("C2",5),rep("C3",5),rep("C4",5),rep("C5",5),rep("C6",5),rep("C7",5),rep("C8",5),rep("C9",5),rep("C10",5),rep("C11",5),rep("C12",5))


# Set the theoretical effects

# Let's have treatment 1 and 2 having no true effect, only treatment 3
effect_treatment_3 = 20
# A fairly large cluster effect
magnitude_cluster_effect=10
# And a fairly small random effect
magnitude_random_effect=3

# Do calculations to obtain the data variable
effect_treatments = rep(0,length(unique(treatment_factor)))
names(effect_treatments)=unique(treatment_factor)
effect_treatments["T3"]=effect_treatment_3

contribution_treatments = effect_treatments[treatment_factor]

effect_cluster=rnorm(length(unique(cluster_factor)),mean=0,sd=magnitude_cluster_effect)
names(effect_cluster)=unique(cluster_factor)

contribution_clusters = effect_cluster[cluster_factor]

contribution_random = rnorm(length(treatment_factor),mean=0,sd=magnitude_random_effect)

data=contribution_random+contribution_clusters+contribution_treatments

names(data)=NULL

get_t_test_matrix_moulton(treatment_factor,data,cluster_factor)
}

\keyword{ misc }

