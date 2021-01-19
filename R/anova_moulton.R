anova_moulton<-function(linmod,data_object,grouping_col,use_predictors=TRUE)
{
    
    mf=model.frame(linmod)
    
    outcome=data_object[,colnames(mf)[1]]
    
    
    
    group = data_object[,grouping_col]
   
   
    anv=anova(linmod)
    
    
    
    moulton_factors = vector(mode="numeric",length=length(rownames(anv))-1)
    
    names(moulton_factors)=rownames(anv)[1:(length(rownames(anv))-1)]
    
    ICC_outcome = ICC(outcome,group)
    
    ICC_predictor = moulton_factors
    
    
    for(moulton_col in names(moulton_factors))
    {
        
        
     
        if(use_predictors) {
         estimator = as.numeric(data_object[,moulton_col]) } else { estimator=NULL }
     
     
        moulton_factors[moulton_col]=moulton_factor(outcome,group,estimator)
        
        ICC_predictor[moulton_col]=ICC(estimator,group)
        
        
    }
    
    
    
    
    
    
    
    anv[1:(dim(anv)[1]-1),"F value"]=anv[1:(dim(anv)[1]-1),"F value"]/moulton_factors^2
    
    # This is the usual number of degrees of freedom. If there is no within cluster variability of the predictors, one should rather use the number of clusters
    
    
    df2=rep(anv["Residuals","Df"],length(moulton_factors))
    
    delta = dim(data_object)[1]-df2[1]
    
    names(df2)=names(moulton_factors)
    
    for(df2_col in names(df2))
    {
        
        predictor = data_object[,df2_col]
        
        OK=!is.na(predictor) & !is.nan(predictor) & !is.na(group) & !is.nan(group)
        
        predictorOK = predictor[OK]
        
        groupOK = group[OK]
        
        
        
        ag=aggregate( predictorOK ~groupOK,FUN=function(x){return(all(duplicated(x)[-1]))})
        
        if(all(ag$predictorOK)) { df2[df2_col]=length(unique(group))-delta  }
    }
    
    
    
    df1=anv[1:(dim(anv)[1]-1),"Df"]
    
    anv[1:(dim(anv)[1]-1),"Pr(>F)"]=pf(anv[1:(dim(anv)[1]-1),"F value"],df1,df2,lower.tail=FALSE)
    
    attr(anv,"moulton")=moulton_factors
    
    attr(anv,"df1")=df1
    
    attr(anv,"df2")=df2
    
    attr(anv,"ICC")=list(outcome=ICC_outcome,predictors=ICC_predictor)
    
    F_overall = mean(anv[1:(dim(anv)[1]-1), "F value"]*df1)/mean(df1)
    
    df1_overall = sum(df1)
    
    df2_overall = mean(df2*df1)/mean(df1)
    
    F_vector_overall = c(F_overall,df1_overall,df2_overall,pf(F_overall,df1_overall,df2_overall,lower.tail=FALSE))
    
    names(F_vector_overall) = c("F value","df1","df2","Pr(>F)")
    
    attr(anv,"overall_fstatistics")=F_vector_overall
    
    # Overall F-statistics, simplified version
    
        
    
    
    return(anv)
    
    
    
}