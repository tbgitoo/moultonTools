ICC<-function(outcome,group,method="unbiased")
{
    # From Mostly Harmless Econometrics, D. Angrist
    
    # First, get the ICC from simple linear regression:
    
    
   
    
    
    
    group=as.factor(group)
    
    OK = !is.na(group) & !is.nan(group) & !is.na(outcome) & !is.nan(outcome)
    
    # Check for intragroup variability, if this is zero
    
    group=group[OK]
    outcome=outcome[OK]
    
    # Test for the trivial case of completely homogeneous groups first:
    
    # For this, all groups should have length larger than 1
    
    group_length=aggregate(outcome ~ group, FUN=length)
    
    if(all(group_length$outcome==1)) { return (1) }
    
    # In R, var(x) gives the sample variance with n-1 degrees of liberty, this is a convenience function for the sum of squares sum(xi-xbar)^2
    sum_squares=function(x)
    {
        return (sum((x-mean(x))^2))
    }
    
    
    # For the unbiased method
    intraclass_sum=function(x,x_mean)
    {
        # The original formula is (xi-x_mean)*(xj-x_mean), i and j being non equal. x_mean is the general means and thus externally imposed (from Mostly harmless econometrics, Angrist and Pischke)
        # One can develop this:
        # Sum xi*xj- n_g*x_mean*Sum xi - ng*x_mean*Sum xj + ng^2 sum x_mean^2 - sum (xi-x_mean)^2 where now the sum is over all combinations possible including i=j
        # Then, sum(xi) becomes ng*xbar where ng is the group size and xbar is the group mean
        # Sum xi*xj - 2*n_g^2*x_mean*x_bar + n_g^2*x_mean^2 - sum(xi^2) + 2*x_mean*sum(xi) - ng*x_mean^2
        # Sum xi*xj - n_g^2*x_mean*(x_bar+x_bar-x_mean)-sum(xi^2) + 2*ng*x_mean*x_bar - ng*x_mean^2
        # Sum xi* sum xj - n_g^2*x_mean*(x_bar+x_bar-x_mean)-sum(xi^2) - n_g*x_mean*(x_mean - x_bar - x_bar)
        # n_g^2*(x_mean-x_bar)^2 - sum(xi-x_bar)^2-ng*(x_mean-x_bar)^2
        # n_g*(n_g-1)*(x_mean-x_bar)^2-sum(xi-x_bar)^2
        # n_g*(n_g-1)*(x_mean-x_bar)^2-var(xi)*(n_g-1) # Sample variance in R with n-1 degrees of freedom
        
        n_g = length(x)
        x_bar = mean(x)
            
        return(n_g*(n_g-1)*(x_mean-x_bar)^2-sum_squares(x))
        
        
        
        
    }
    
    # The outcome is numerical, apply the standard formulae
    
    if(!is.factor(outcome))
    {
    
    x_mean = mean(outcome)
    
    if(method=="unbiased" | method=="Fisher") #from Mostly harmless econometrics, Angrist and Pischke, eq. 8.2.5, transformed for more compact evaluation
    {
    
        ng=aggregate(outcome~group,FUN=length)$outcome
        zg=aggregate(outcome~group,FUN=mean)$outcome
        n=length(outcome)
        return((sum(ng^2*(zg-mean(outcome))^2)/(sum((outcome-mean(outcome))^2)/n)-n)/sum(ng*(ng-1)))
    
    } else if(method=="unbiased_unequal_cluster_size") {
        # We find the formula in Mostly harmless econometric to not give exactly 1 for homogeneous cluster (no intra-cluster variability) when the clusters
        # are not of equal size. To address this problem, the formula here can be used. It is engineered to give approximatly 0 for non-clustered variables,
        # and exactly 1 for perfect clustering
        ng=aggregate(outcome~group,FUN=length)$outcome
        N=length(ng)
        zg=aggregate(outcome~group,FUN=mean)$outcome
        n=length(outcome)
        
        return((n-3)/n/(n-N-2)*(sum(ng*(zg-mean(outcome))^2)/(sum((outcome-mean(outcome))^2)/n)-n*(N-1)/(n-3)))
        
        #return((sum(ng*(zg-mean(z))^2)/(sum((z-mean(z))^2)/n)))
        
    
    } else if(method=="ANOVA")
    { # Group variance over total residual variance, general method, for example eq. 8.2.3 in Mostly harmless econometrics, Angrist and Pischke
    
    ag=aggregate( outcome ~group,FUN=function(x){return(all(duplicated(x)[-1]))})
    
    
    # Are all groups homogeneous
    if(all(ag$outcome)) {return(1)}
    
    outcome=as.numeric(outcome)
    
    
    linmod = lm(outcome~group)
    
    Var_intragroup = anova(linmod)["Residuals","Sum Sq"]
    Var_intergroup = anova(linmod)["group","Sum Sq"]
    
    ICC_m=Var_intergroup/(Var_intergroup+Var_intragroup)
    
    return(ICC_m)
    
    } else {
        
        stop(paste("Method ", method," not implemented", sep=""))
    }
    
    

    } else { # If the outcome is a factor, we assume this is non-ordered. We decompose this to 0-1 variables for each level and sum up
        
        sum_ICC=0
        n_ICC=0
        for(theLevel in levels(outcome))
        {
            levelOutcome = as.numeric(outcome==theLevel)
            sum_ICC=sum_ICC+ICC(levelOutcome,group,method=method)
            n_ICC=n_ICC+1
        }
    
        return (sum_ICC/n_ICC)
    
        
    }
    
	
}

