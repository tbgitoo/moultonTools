moulton_factor<-function(outcome,group,estimator=NULL,...)
{
    # From Mostly Harmless Econometrics: An Empiricist's Compagnion, Joshua D. Angrist and Joern-Steffen Pischke
    
    # First, get the ICC from simple linear regression:
    
    
    
    
    group=as.factor(group)
    
    intraclass=0
    
    arguments=list(...)
    if("method" %in% names(arguments))
    {
        if(arguments[["method"]]=="Angrist_mixed")
        {
            call_args=c(list("outcome"=outcome,"group"=group))
            for(theArg in names(arguments))
            {
                if(theArg!="method")
                {
                call_args[[theArg]]=arguments[[theArg]]
                } else {
                    call_args[["method"]]="ANOVA"
                }
            }
            
            
            intraclass=do.call(ICC,call_args)
        } else {
            intraclass = ICC(outcome,group,...)
        }
    } else {
        intraclass = ICC(outcome,group,...)
    }
    
    
    
    
    
    n=aggregate(outcome~group,FUN=length)
    
    px = 1
    
    if(is.null(estimator)) # If no estimator is provided, we assume no intragroup  variability
    {
        
        # Groups of size 1 do not contribute to the Moulton factor
        
       px=1
        
    } else {
       # Allow here for Angrist's mixed method: Fisher method for px but variance method for the ICC
       arguments=list(...)
       if("method" %in% names(arguments))
       {
           if(arguments[["method"]]=="Angrist_mixed")
           {
               call_args=c(list("outcome"=estimator,"group"=group),arguments)
               call_args[["method"]]="unbiased"
               call_args[["var_n"]]="n-1"
               
               px=do.call(ICC,call_args)
           } else {
               px=ICC(estimator,group,...)
           }
       } else {
           px=ICC(estimator,group,...)
       }
       
        
        
    }
    
    
    
    var_n_outcome = sum((n$outcome-mean(n$outcome))^2)/length(n$outcome)
    
    mf=(1+(var_n_outcome/mean(n$outcome)+mean(n$outcome)-1)*px*intraclass)
    
    # If all outcomes are identical, in principle, we cannot estimate variance. Return the default value then
    if(all(outcome==outcome[1])) {mf=1}
    
    if(mf<1) {mf=1}
    
    ret=sqrt(mf)
    
    attr(ret,"ICC")=intraclass
    
    attr(ret,"px")=px
    
    return(ret)
    
    

    
    
	
}
