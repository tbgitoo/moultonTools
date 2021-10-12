df_t_cluster<-function (x, y = NULL,cluster_x=1:length(x),cluster_y=1:length(y),var.equal=TRUE,...)
{
    
    
    if(is.null(y))
    {
        ICC_val=ICC(x,cluster_x,...)
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        df_clustered=length(unique(cluster_x))-1
        df_unclustered=length(x)-1
        df=(1/df_unclustered+ICC_val^2*(1/df_clustered-1/df_unclustered))^(-1)
        return(df)
    } else {
        ICC_val=ICC(c(x-mean(x),y-mean(y)),c(cluster_x,cluster_y),...)
        if(ICC_val <0 ) {ICC_val = 0}
        if(ICC_val > 1) {ICC_val = 1}
        if(var.equal)
        {
          df_clustered=length(unique(cluster_x))+length(unique(cluster_y))-2
          df_unclustered=length(x)+length(y)-2
        } else {
            stderrx <- sqrt(var(x)/length(x))
            stderry <- sqrt(var(y)/length(y))
            stderr=sqrt(stderrx^2+stderry^2)
            cluster_nx=length(unique(cluster_x))
            cluster_ny=length(unique(cluster_y))
            
            df_clustered <- stderr^4/(stderrx^4/(cluster_nx - 1) + stderry^4/(cluster_ny - 1))
            df_unclustered <- stderr^4/(stderrx^4/(length(x) - 1) + stderry^4/(length(y) - 1))
            
        }
        df=(1/df_unclustered+ICC_val^2*(1/df_clustered-1/df_unclustered))^(-1)
        return(df)
        
    }
}
