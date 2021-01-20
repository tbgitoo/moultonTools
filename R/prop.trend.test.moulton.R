prop.trend.test.moulton<-function (binaryResponse, score, cluster=seq_along(binaryResponse),... )
{
    
    x <- as.vector(table(as.logical(binaryResponse),as.numeric(score))["TRUE",])
    n <- as.vector(table(as.numeric(score)))
    
    scores <- as.vector(as.numeric(dimnames(table(as.numeric(score)))[[1]]))
    
    test_result=prop.trend.test(x,n,scores)
    
    mf <- moulton_factor(binaryResponse,cluster,score,...)
    
    test_result[["statistic"]][["X-squared"]]=test_result[["statistic"]][["X-squared"]]/mf^2
    test_result[["p.value"]]=pchisq(as.numeric( test_result[["statistic"]][["X-squared"]]), 1, lower.tail = FALSE)
    test_result[["method"]] = paste(test_result[["method"]],", with Moulton correction",sep="")
    test_result[["moulton_factor"]] <- mf
    
    return(test_result)
}
