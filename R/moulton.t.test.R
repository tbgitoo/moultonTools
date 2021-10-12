moulton.t.test<-function (x, y = NULL,cluster_x=1:length(x),cluster_y=NULL,method_df="cluster", alternative = c("two.sided", "less", "greater"),
mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
...)
{
    # First, if there are no clustering variables provided, this is assumed to be a standard t-test
    
    if(is.null(cluster_x) & is.null(cluster_y)) { return(t.test(x, y, alternative = alternative,
        mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level,
        ...))}
    
    if(is.null(y) & is.null(cluster_x)) { return(t.test(x, y, alternative = alternative,
        mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level,
        ...))}
    
    
    # If clustering information is provided for only one of the two variables, the other is supposed to unclustered
    
    if(is.null(cluster_x)) { cluster_x = 1:length(x) }
    if(is.null(cluster_y) & !is.null(y)) {cluster_y = length(cluster_x)+1:length(y)}
    
    
    
    
    # Do t-test standard stuff from R standard implementation
    
    alternative <- match.arg(alternative)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
    conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        if (paired)
        xok <- yok <- complete.cases(x, y)
        else {
            yok <- !is.na(y)
            xok <- !is.na(x)
        }
        y <- y[yok]
    }
    else {
        dname <- deparse(substitute(x))
        if (paired)
        stop("'y' is missing for paired test")
        xok <- !is.na(x)
        yok <- NULL
    }
    x <- x[xok]
    if (paired) {
        x <- x - y
        y <- NULL
    }
    nx <- length(x)
    cluster_nx <- length(unique(cluster_x))
    mx <- mean(x)
    vx <- var(x)
    mf <- NA
    if (is.null(y)) {
        if (nx < 2)
        stop("not enough 'x' observations")
        
        
        
        if(cluster_nx==1)
        {
        
            stop("All observations in a single cluster, cannot perform testing")
            
        }
        
        # Conservative estimate of the degree of freedom: the number of clusters
        
        df <- cluster_nx - 1
        
        if(method_df=="ICC")
        {
            df=df_t_cluster(x=x,cluster_x=cluster_x,...)
            
        }
        
        moulton_factor_x=1
        
        arguments=list(...) # Potentially, the users wants a specific ICC method
        if("method" %in% names(arguments))
        {
        
        moulton_factor_x = moulton_factor(x,cluster_x,method=arguments["method"])
        
        } else {
            moulton_factor_x=moulton_factor(x,cluster_x)
        }
        
        
        # Correct the variance by multiplying with the Moulton factor squared
        
        vx<-vx*moulton_factor_x^2
        mf<-moulton_factor_x
        
        stderr <- sqrt(vx/nx)
        
        if (stderr < 10 * .Machine$double.eps * abs(mx))
        stop("data are essentially constant")
        
        
        tstat <- (mx - mu)/stderr
        
        
        
        method <- if (paired)
        "Paired t-test"
        else "One Sample t-test, with Moulton correction"
        estimate <- setNames(mx, if (paired)
        "mean of the differences, with Moulton correction"
        else "mean of x, with Moulton correction")
    }
    else {
        ny <- length(y)
        cluster_ny <- length(unique(cluster_y))
        if (nx < 1 || (!var.equal && nx < 2))
        stop("not enough 'x' observations")
        if (ny < 1 || (!var.equal && ny < 2))
        stop("not enough 'y' observations")
        if (var.equal && nx + ny < 3)
        stop("not enough observations")
        
        if(var.equal && cluster_ny+cluster_nx < 3)
        stop("Too few clusters in x and y")
        
        my <- mean(y)
        vy <- var(y)
        method <- paste(if (!var.equal)
        "Welch", "Two Sample t-test")
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        
        mf<-1
        
        arguments=list(...) # Potentially, the users wants a specific ICC method
        if("method" %in% names(arguments))
        {
        
        mf<-moulton_factor(c(x,y),c(cluster_x,cluster_y),c(rep(1,nx),rep(2,ny)),
            method=arguments["method"])
        
        } else {
            mf<-moulton_factor(c(x,y),c(cluster_x,cluster_y),c(rep(1,nx),rep(2,ny)))
        }
        
        
        if (var.equal) {
            df <- nx + ny - 2
            v <- 0
            if (nx > 1)
            v <- v + (nx - 1) * vx
            if (ny > 1)
            v <- v + (ny - 1) * vy
            v <- v/df
            
            
            v<-v*mf
            
            
            
            stderr <- sqrt(v * (1/nx + 1/ny))
            
            # Correct for clustering in the degrees of freedom
            
            df<-cluster_nx + cluster_ny -2
            
            
            if(method_df=="ICC")
            {
                df=df_t_cluster(x=x,cluster_x=cluster_x,y=y,cluster_y=cluster_y,var.equal=TRUE,...)
                
            }
            
            
        } else {
            
            stderrx <- sqrt(vx/nx)
            stderry <- sqrt(vy/ny)
            
            
           
            
            # df <- cluster_nx + cluster_ny -2
            
             stderr <- sqrt(stderrx^2 + stderry^2)
            df <- stderr^4/(stderrx^4/(cluster_nx - 1) + stderry^4/(cluster_ny -
            1))
             if(method_df=="ICC")
             {
                 df=df_t_cluster(x=x,cluster_x=cluster_x,y=y,cluster_y=cluster_y,var.equal=FALSE,...)
             }
             stderr <- stderr*mf
        }
        if (stderr < 10 * .Machine$double.eps * max(abs(mx),
        abs(my)))
        stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
        pval <- pt(tstat, df)
        cint <- c(-Inf, tstat + qt(conf.level, df))
    }
    else if (alternative == "greater") {
        pval <- pt(tstat, df, lower.tail = FALSE)
        cint <- c(tstat - qt(conf.level, df), Inf)
    }
    else {
        pval <- 2 * pt(-abs(tstat), df)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y))
    "difference in means"
    else "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
    conf.int = cint, estimate = estimate, null.value = mu,
    alternative = alternative, method = method, data.name = dname,moulton_factor=mf)
    class(rval) <- "htest"
    return(rval)
}
