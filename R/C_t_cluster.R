C_t_cluster<-function (x, y = NULL,cluster_x=1:length(x),cluster_y=1:length(y),var.equal=TRUE,method="ANOVA",...)
{
    
    restrict_ICC=FALSE
    if(is.null(y))
    {
        ICC_val=ICC(x,cluster_x,method=method,...)
        if(restrict_ICC)
        {
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        }
        C_clustered=(length(x)-max(table(cluster_x)))/length(x)
        C_unclustered=(length(x)-1)/length(x)
        
    } else {
        ICC_val=ICC(c(x-mean(x),y-mean(y)),c(cluster_x,cluster_y),method=method,...)
        if(restrict_ICC)
        {
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        }
        if(var.equal) # This is the homoscedastic case
        {
          x_cluster_size = as.vector(table(cluster_x))
          y_cluster_size = as.vector(table(cluster_y))
          n=sum(x_cluster_size)+sum(y_cluster_size)
          # Take into account unequal sizes of the clusters
          G=length(unique(c(cluster_x,cluster_y))) # Total number of clusters
          # - the two largest clusters for a worst case scenario
          C_clustered=(length(x)+length(y))/(length(x)+length(y)-sum(sort(table(c(cluster_x,cluster_y)),decreasing=TRUE)[1:2]))
          C_unclustered=(length(x)+length(y))/(length(x)+length(y)-2)
        } else {
            # Unequal variances, use the full expression by estimating the variance matrix
            # as contributions from blockwise correlated and uncorrelated variables
            # Regressor weights proportional to variance contribution and inversely proportional to the n's
            X<-c(rep(sqrt(var(x)/length(x)),length(x)),-rep(sqrt(var(y)/length(y)),length(y)))
            Vz=diag(length(X))
            
            Vy=Vz
            for(theCluster in unique(c(cluster_x,cluster_y)))
            {
                Vy[theCluster==c(cluster_x,cluster_y),theCluster==c(cluster_x,cluster_y)]=1
            }
            
            Vz[1:length(x),]=Vz[1:length(x),]*sqrt(var(x))
            Vz[,1:length(x)]=Vz[,1:length(x)]*sqrt(var(x))
            Vz[length(x)+(1:length(y)),]=Vz[length(x)+(1:length(y)),]*sqrt(var(y))
            Vz[,length(x)+(1:length(y))]=Vz[,length(x)+(1:length(y))]*sqrt(var(y))
            Vy[1:length(x),]=Vy[1:length(x),]*sqrt(var(x))
            Vy[,1:length(x)]=Vy[,1:length(x)]*sqrt(var(x))
            Vy[length(x)+(1:length(y)),]=Vy[length(x)+(1:length(y)),]*sqrt(var(y))
            Vy[,length(x)+(1:length(y))]=Vy[,length(x)+(1:length(y))]*sqrt(var(y))
            
            V=(1-ICC_val)*Vz+ICC_val*Vy # Estimate of the variance matrix
            
            # Projection matrix
            
            M=diag(length(X))
            M[1:length(x),1:length(x)]=M[1:length(x),1:length(x)]-1/length(x)
            M[length(x)+(1:length(y)),length(x)+(1:length(y))]=M[length(x)+(1:length(y)),length(x)+(1:length(y))]-1/length(y)
            
            C=(t(X)%*%diag(diag(V))%*%X)/(t(X)%*%diag(diag(M%*%V))%*%X)
            
            return(C)
            
             
            
            
        }
        C=(1/C_unclustered+ICC_val*(1/C_clustered-1/C_unclustered))^(-1)
        return(C)
        
    }
}
