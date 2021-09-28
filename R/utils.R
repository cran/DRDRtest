rrademachar <- function(n){
    ## intput:
    ## n: number of observations
    ##
    ## output:
    ## a vector containing the simulated observations        
    x <- c(1+sqrt(5),1-sqrt(5))/2
    p <- (sqrt(5)-1)/(2*sqrt(5))
    return(sample(x,size=n,replace = TRUE,prob=c(p,1-p)))
}

rtwopoint <- function(n){
    ## intput:
    ## n: number of observations
    ##
    ## output:
    ## a vector containing the simulated observations        
    x <- c(1+sqrt(5),1-sqrt(5))/2
    p <- (sqrt(5)-1)/(2*sqrt(5))
    return(sample(x,size=n,replace = TRUE,prob=c(p,1-p)))
}

wildboot <- function(x,y,h,b,xrange,dist='TwoPoint',x.grid.size){
    n <- length(y)
    yhat <- mean(y)
    res <- y-yhat
    if(dist=='TwoPoint'){
        estar <- matrix(rtwopoint(length(y)*b),nrow=n)*res
    }else if(dist=='Rademachar'){
        estar <- matrix(rrademachar(length(y)*b),nrow=n)*res
    }else{
        stop("dist should be either 'TwoPoint' or 'Rademachar'")
    }
    yboot <- yhat+estar
    tboot.list <- rep(0,length.out=b)
    for(i in 1:b){
        loc.wild <- KernSmooth::locpoly(x,yboot[,i],drv=0,bandwidth = h, range.x=xrange,gridsize = x.grid.size)
        linear.est <- mean(yboot[,i])
        Tboot <- sum((loc.wild$y-linear.est)^2)*(loc.wild$x[2]-loc.wild$x[1])*n*sqrt(h)
        tboot.list[i] <- Tboot
    }
    loc.obs <- KernSmooth::locpoly(x,y,drv=0,bandwidth =h, range.x=xrange,gridsize = x.grid.size)
    Tobs <- sum((loc.obs$y-yhat)^2)*(loc.obs$x[2]-loc.obs$x[1])*n*sqrt(h)

    return(list(T.obs = Tobs, T.boot = tboot.list,loc.fit = loc.obs))
}


wildboot_conditional <- function(xlist,ylist,h,b,xrange, dist='TwoPoint', x.grid.size=401){
    loc.all <- KernSmooth::locpoly(unlist(xlist),unlist(ylist),bandwidth=h, range.x=xrange, gridsize=x.grid.size)
    yboot_list <- list()
    for(c in 1:length((xlist))){
        res <- ylist[[c]]-stats::approx(loc.all,xout=xlist[[c]])$y
        if(dist=='TwoPoint'){
            estar <- matrix(rtwopoint(length(res)*b),nrow=length(res))*res
        }else if(dist=='Rademachar'){
            estar <- matrix(rrademachar(length(res)*b),nrow=length(res))*res
        }else{
            stop("dist should be either 'TwoPoint' or 'Rademachar'")
        }
        yboot_list[[c]] <- estar+stats::approx(loc.all,xout=xlist[[c]])$y
    }

    fit_at_grid <- matrix(0,ncol=length(xlist),nrow=length(loc.all$x))
    ## observed test statistc
    Tobs <- 0
    for(c in 1:length(xlist)){
        loc.temp <- KernSmooth::locpoly(xlist[[c]],ylist[[c]],bandwidth=h,range.x=xrange, gridsize = x.grid.size)
        fit_at_grid[,c] <- loc.temp$y
        if(c>=2){
            Tobs <- Tobs +sum((fit_at_grid[,1:(c-1)]-loc.temp$y)^2)*(loc.temp$x[2]-loc.temp$x[1])                
        }
    }
    
    ## Bootstrap test statistic
    tboot.list <- rep(0, length.out=b)
    for(i in 1:b){
        Tboot<-0
        for(c in 1:length(xlist)){
            loc.temp <- KernSmooth::locpoly(xlist[[c]],yboot_list[[c]][,i],bandwidth=h,range.x=xrange, gridsize = x.grid.size)
            fit_at_grid[,c] <- loc.temp$y
            if(c>=2){
                Tboot <- Tboot +sum((fit_at_grid[,1:(c-1)]-loc.temp$y)^2)*(loc.temp$x[2]-loc.temp$x[1])                
                }
        }
        tboot.list[i] <- Tboot
    }
    return(list(T.obs = Tobs, T.boot = tboot.list))
}
