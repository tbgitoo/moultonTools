df_t_cluster<-function (x, y = NULL,cluster_x=1:length(x),cluster_y=1:length(y),var.equal=TRUE,method="ANOVA",...)
{
    
    
    if(is.null(y))
    {
        ICC_val=ICC(x,cluster_x,method=method,...)
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        df_clustered=length(unique(cluster_x))-1
        df_unclustered=length(x)-1
        df=(1/df_unclustered+ICC_val^2*(1/df_clustered-1/df_unclustered))^(-1)
        return(df)
    } else {
        ICC_val=ICC(c(x-mean(x),y-mean(y)),c(cluster_x,cluster_y),method=method,...)
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        if(var.equal) # This is the homoscedastic case
        {
          x_cluster_size = as.vector(table(cluster_x))
          y_cluster_size = as.vector(table(cluster_y))
          n=sum(x_cluster_size)+sum(y_cluster_size)
          # Take into account unequal sizes of the clusters
          G=length(unique(c(cluster_x,cluster_y))) # Total number of clusters
          df_clustered=n^2/(sum(x_cluster_size^2)+sum(y_cluster_size^2))*(G-2)/G
          df_unclustered=length(x)+length(y)-2
        } else {
            stderrx <- sqrt(var(x)/length(x))
            stderry <- sqrt(var(y)/length(y))
            stderr=sqrt(stderrx^2+stderry^2)
            df_unclustered <- stderr^4/(stderrx^4/(length(x) - 1) + stderry^4/(length(y) - 1))
            # Clustered case with unequal variance
            x_cluster_size = as.vector(table(cluster_x))
            y_cluster_size = as.vector(table(cluster_y))
            cluster_nx=length(unique(cluster_x))
            cluster_ny=length(unique(cluster_y))
            df_x = sum(x_cluster_size)^2/(sum(x_cluster_size^2))-1
            df_y = sum(y_cluster_size)^2/(sum(y_cluster_size^2))-1
            # Take into account unequal sizes of the clusters
            
            df_clustered <- stderr^4/(stderrx^4/df_x + stderry^4/df_y)
            
            
            
            
        }
        df=(1/df_unclustered+ICC_val^2*(1/df_clustered-1/df_unclustered))^(-1)
        return(df)
        
    }
}
