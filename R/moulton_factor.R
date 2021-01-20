moulton_factor<-function(outcome,group,estimator=NULL,...)
{
    # From Mostly Harmless Econometrics: An Empiricist's Compagnion, Joshua D. Angrist and Joern-Steffen Pischke
    
    # First, get the ICC from simple linear regression:
    
    
    
    
    group=as.factor(group)
    
    intraclass = ICC(outcome,group,...)
    
    
    
    n=aggregate(outcome~group,FUN=length)
    
    px = 1
    
    if(is.null(estimator)) # If no estimator is provided, we assume no intragroup  variability
    {
        
        # Groups of size 1 do not contribute to the Moulton factor
        
       px=1
        
    } else {
        
        
       px=ICC(estimator,group,...)
        
       
        
        
    }
    
    var_n_outcome = sum((n$outcome-mean(n$outcome))^2)/length(n$outcome)
    
    mf=(1+(var_n_outcome/mean(n$outcome)+mean(n$outcome)-1)*px*intraclass)
    if(mf<1) {mf=1}
    
    ret=sqrt(mf)
    
    attr(ret,"ICC")=intraclass
    
    attr(ret,"px")=px
    
    return(ret)
    
    

    
    
	
}
