#' The base function for testing effect modifiers
#'
#' This is the base function for testing whether a discrete covariate is an effect modifier.
#' 
#' @param ylist A list containing vectors of outcomes for each class
#' @param alist A list containing vectors of  treatment levels (dosage) for each class
#' @param pilist A list containing vectors of propensity scores for each class
#' @param varpilist A list containing vectors of mean propensity scores for each class
#' @param mulist A list containing vectors of outcome regression function values for each class
#' @param malist A list containing vectors of  mean outcome regression values for each class
#' @param arange A vector of length 2 giving the lower bound and upper bound of  treatment levels
#' @param h bandwidth to be used in kernel regression. If not specified, will by default use "rule of thumb" bandwidth selector
#' @param b number of Bootstrap samples to be generated
#' @param dist distibution used to generate residuals for Bootstrap samples. Currently only have two options, "TwoPoint" and "Rademachar"
#' @param a.grid.size size of equally spaced grid points over \code{arange} to be generate for numerically evaluating the integral in test statistic
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}     
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' 
#' @export
#' @examples
#' d <- 4
#' n <- 200
#' sigma <- 0.5
#' delta <- 1
#' height <-1
#' arange <- c(0,5)
#' triangle <- function(a,height){
#'    y <- exp(-a^2/((1/2)^2))*height
#'    return(y)
#'}
#'mu.mod<-function(a,l,delta,height){
#'    mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1*delta))+
#'          triangle(a-2.5,height)+a*(-0.1*l[,1]+0.1*delta*l[,4])
#'    return(mu)
#'}
#'l <- matrix(rnorm(n*d),ncol=d)
#'l[,4] <- ifelse(l[,4]>0,1,0)
#'colnames(l) <- paste("l",1:4,sep="")
#'
#'logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0))
#'lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#'a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#'mu <- mu.mod(a,l,delta,height)
#'residual.list <- rnorm(n,mean=0,sd =sigma)
#'y <- mu+residual.list
#'
#'class_label <- l[,4]
#' ylist <- split(y,class_label)
#'alist <- split(a,class_label)
#'pilist <- split(pmin(dbeta(a/5,shape1=lambda,shape2=1-lambda)/5,100),class_label)
#'mulist <- split(mu,class_label)
#'
#'varpilist <- list()
#'malist <- list()
#'for(c in c(0,1)){
#'    ac <- a[class_label==c]
#'    lc <- l[class_label==c,]
#'
#'    logit.lambdac <- as.numeric(lc[rep(1:nrow(lc),nrow(lc)),]%*%c(0.1,0.1,-0.1,0))
#'    lambdac <- exp(logit.lambdac)/(1+exp(logit.lambdac))
#'    varpic <- colMeans(matrix(pmin(dbeta(rep(ac,each=length(ac))/5,
#'                                   shape1=lambdac,
#'                                   shape2 = 1-lambdac)/5,100),nrow=length(ac)))
#'
#'    mac <- colMeans(matrix(mu.mod(rep(ac,each=length(ac)),
#'                                      lc[rep(1:nrow(lc),nrow(lc)),],
#'                                      delta,height),
#'                           nrow=length(ac)))
#'
#'    varpilist[[as.character(c)]]<-varpic
#'    malist[[as.character(c)]] <- mac
#'    }
#'    
#'out <- drdrtest_em.base(ylist,alist,pilist,varpilist,mulist,malist,arange)

drdrtest_em.base <-function(ylist,alist,pilist,varpilist,mulist,malist,arange,h=NULL,b=1000,dist='TwoPoint',a.grid.size=401){
    pseudo_list <- list()
    for(i in 1:length(ylist)){
        pseudo.out <- (ylist[[i]]-mulist[[i]])/(pilist[[i]]/varpilist[[i]])+malist[[i]]
        pseudo_list[[i]] <- pseudo.out
        }
    n <- length(unlist(ylist))
    if(is.null(h)){
        h <- 0.9*n^(-1/5)*min(stats::sd(unlist(alist)),stats::IQR(unlist(alist))/1.34)
    }
    test.out <- wildboot_conditional(alist,pseudo_list,h,b,arange,dist,x.grid.size=a.grid.size)
    return(list(p.value = mean(test.out$T.boot>=test.out$T.obs),
                test.stat = test.out$T.obs,
                Bootstrap.samples = test.out$T.boot,
                bandwidth=h))
}
#' The base function for testing a effect modifier with user specified nuisance functions
#'
#' This is the  function for testing whether a discrete covariate is an effect modifier with user specified nuisance functions
#' 
#' @param y A vector containing the outcomes for each observation
#' @param a A vector containing the treatment levels (dosage) for each observation
#' @param l A data.frame containing the observations of covariates
#' @param class_label A vector containing the class label (label for the effect modifier) for each observation.
#' @param arange A vector of length 2 giving the lower bound and upper bound of  treatment levels
#' @param pifunc A user specifid function or wapper that takes treatment a as the first argument and
#' covariates l as the second argument and return propensit scores
#' @param mufunc A user specifid function or wapper that takes treatment a as the first argument and
#' covariates l as the second argument and return outcome regression values
#' @param h bandwidth to be used in kernel regression. If not specified, will by default use "rule of thumb" bandwidth selector
#' @param b number of Bootstrap samples to be generated
#' @param dist distibution used to generate residuals for Bootstrap samples. Currently only have two options, "TwoPoint" and "Rademachar"
#' @param pi.low Lower bound to truncate propensity scores
#' @param a.grid.size size of equally spaced grid points over \code{arange} to be generate for numerically evaluating the integral in test statistic
#' 
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}     
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' @export
#' @examples
#' d <- 4
#' n <- 200
#' sigma <- 0.5
#' delta <- 1
#' height <-1
#' arange <- c(0,5)
#' triangle <- function(a,height){
#'    y <- exp(-a^2/((1/2)^2))*height
#'    return(y)
#'}
#'mu.mod<-function(a,l,delta,height){
#'    mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1*delta))+
#'          triangle(a-2.5,height)+a*(-0.1*l[,1]+0.1*delta*l[,4])
#'    return(mu)
#'}
#'l <- matrix(rnorm(n*d),ncol=d)
#'l[,4] <- ifelse(l[,4]>0,1,0)
#'colnames(l) <- paste("l",1:4,sep="")
#'
#'logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0))
#'lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#'a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#'mu <- mu.mod(a,l,delta,height)
#'residual.list <- rnorm(n,mean=0,sd =sigma)
#'y <- mu+residual.list
#'
#'class_label <- l[,4]
#'
#' pifunc <- function(a,l){
#'    l <- as.matrix(l)
#'    logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0))
#'    lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#'    return(pmin(dbeta(a/5,shape=lambda,shape2=1-lambda)/5,100))    
#'}
#'
#'mufunc <- function(a,l){
#'    return(mu.mod(a,as.matrix(l),delta,height))
#'}
#'
#'out <- drdrtest_em(y,a,l,class_label,arange,pifunc,mufunc)

drdrtest_em <- function(y,a,l,class_label,arange,pifunc,mufunc,h=NULL,b=1000,dist='TwoPoint',pi.low=0.01,a.grid.size=401){
    approx.fn <- function(x,y,z){
        stats::predict(stats::smooth.spline(x,y),x=z)$y
    }

    class_label.factor <- factor(class_label,levels = sort(unique(class_label)))
    la <- data.frame(l,a=a)
    a.vals <-seq(arange[1],arange[2],length.out=a.grid.size)
    
    pihat <- pmax(pifunc(a,l),pi.low)
    pihat_list <- split(pihat, class_label.factor)

    muhat <- mufunc(a,l)
    muhat_list <- split(muhat, class_label.factor)

    varpihat_list <-list()
    mahat_list<- list()
    for(c in sort(unique(class_label))){
        lc <- l[class_label==c,]
        ac <- a[class_label==c]
        la.new <- data.frame(lc[rep(1:nrow(lc),length(a.vals)),],a=rep(a.vals,each=nrow(lc)))
        pihat.mat <- matrix(pmax(pifunc(la.new[,ncol(la.new)],la.new[,-ncol(la.new)]),pi.low),ncol=length(a.vals))
        varpihat <- approx.fn(a.vals,apply(pihat.mat,2,mean),ac)
        varpihat_list[[as.character(c)]] <- varpihat

        muhat.mat <- matrix(mufunc(la.new[,ncol(la.new)],la.new[,-ncol(la.new)]),ncol=length(a.vals))
        mahat_list[[as.character(c)]] <- approx.fn(a.vals,apply(muhat.mat,2,mean),ac)
    }
    ylist <- split(y,class_label.factor)
    alist <- split(a,class_label.factor)
    return(drdrtest_em.base(ylist,alist,pihat_list,varpihat_list,
                            muhat_list,mahat_list,arange, h,b,dist, a.grid.size))
    
}

#' The function for testing a effect modifier with SuperLearner
#'
#' This is the  function for testing whether a discrete covariate is an effect modifier with SuperLearner
#' 
#' @param y A vector containing the outcomes for each observation
#' @param a A vector containing the treatment levels (dosage) for each observation
#' @param l A data.frame containing the observations of covariates
#' @param class_label A vector containing the class label (label for the effect modifier) for each observation.
#' @param arange A vector of length 2 giving the lower bound and upper bound of  treatment levels
#' @param pi.sl.lib Models will be used by SuperLearner to estiamte propensity scores
#' @param mu.sl.lib Models will be used by SuperLearner to estiamte outcome regression function
#' @param mu.family Type of response. Currently only support "gaussian" and "binomial"
#' @param h bandwidth to be used in kernel regression. If not specified, will by default use "rule of thumb" bandwidth selector
#' @param b number of Bootstrap samples to be generated
#' @param dist distibution used to generate residuals for Bootstrap samples. Currently only have two options, "TwoPoint" and "Rademachar"
#' @param pi.low Lower bound to truncate propensity scores
#' @param pi.var.low Lower bound to truncate conditional variance of treament (used in propensity score estimation).
#' @param a.grid.size size of equally spaced grid points over \code{arange} to be generate for numerically evaluating the integral in test statistic
#' 
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}     
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' @export
#' @examples
#' d <- 4
#' n <- 200
#' sigma <- 0.5
#' delta <- 1
#' height <-1
#' arange <- c(0,5)
#' triangle <- function(a,height){
#'    y <- exp(-a^2/((1/2)^2))*height
#'    return(y)
#'}
#'mu.mod<-function(a,l,delta,height){
#'    mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1*delta))+
#'          triangle(a-2.5,height)+a*(-0.1*l[,1]+0.1*delta*l[,4])
#'    return(mu)
#'}
#'l <- matrix(rnorm(n*d),ncol=d)
#'l[,4] <- ifelse(l[,4]>0,1,0)
#'colnames(l) <- paste("l",1:4,sep="")
#'
#'logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0))
#'lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#'a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#'mu <- mu.mod(a,l,delta,height)
#'residual.list <- rnorm(n,mean=0,sd =sigma)
#'y <- mu+residual.list
#'
#'class_label <- l[,4]
#'out <- drdrtest_em.superlearner(y,a,l,l[,4],arange,pi.sl.lib=c("SL.glm"),mu.sl.lib=c("SL.glm"))
drdrtest_em.superlearner <- function(y,a,l,class_label, arange,pi.sl.lib=c("SL.earth","SL.glm","SL.gam","SL.glmnet"),
    mu.sl.lib=c("SL.earth","SL.glm","SL.gam","SL.glmnet"),
    mu.family ="gaussian",h=NULL,b=1000,
    dist='TwoPoint',pi.low=0.01,pi.var.low=0.01,a.grid.size=401){
    approx.fn <- function(x,y,z){
        stats::predict(stats::smooth.spline(x,y),x=z)$y
    }

    class_label.factor <- factor(class_label,levels = sort(unique(class_label)))
    la <- data.frame(l,a=a)
    a.vals <-seq(arange[1],arange[2],length.out=a.grid.size)

    ## calculating pihat
    pimod <- SuperLearner::SuperLearner(Y=a,X=as.data.frame(l),SL.library=pi.sl.lib,newX=la[,-ncol(la)])
    pimod.fit <- as.numeric(pimod$SL.predict)
    
    sq.res <- (a-pimod.fit)^2
    pimod2 <- SuperLearner::SuperLearner(Y=sq.res,X=la[,-ncol(la)],SL.library=pi.sl.lib,newX=la[,-ncol(la)])
    pimod2.fit <- pmax(as.numeric(pimod2$SL.predict), pi.var.low)

    a.std <- (a-pimod.fit)/sqrt(pimod2.fit)
    pihat <- pmax(approx.fn(stats::density(a.std)$x, stats::density(a.std)$y, a.std)/sqrt(pimod2.fit),pi.low)
    pihat_list <- split(pihat,class_label.factor)

    ## calcualte muhat
    mumod <- SuperLearner::SuperLearner(Y=y,X=data.frame(l,a=a),SL.library=mu.sl.lib,newX=la,family=mu.family)
    muhat <- as.numeric(mumod$SL.predict)
    muhat_list <- split(muhat,class_label.factor)
    
    ## calculating varpihat
    varpihat_list <-list()
    mahat_list<- list()
    for(c in sort(unique(class_label))){
        lc <- l[class_label==c,]
        ac <- a[class_label==c]
        la.new <- data.frame(lc[rep(1:nrow(lc),length(a.vals)),],a=rep(a.vals,each=nrow(lc)))
        pimod.fit.anew <- as.numeric(stats::predict(pimod,newdata=la.new[,-ncol(la.new)])$pred)
        pimod2.fit.anew <- pmax(as.numeric(stats::predict(pimod2,newdata=la.new[,-ncol(la.new)])$pred),pi.var.low)
        anew.std <- (la.new$a-pimod.fit.anew)/sqrt(pimod2.fit.anew)
        pihat.mat <- matrix(approx.fn(stats::density(a.std)$x, stats::density(a.std)$y, anew.std)/sqrt(pimod2.fit.anew),ncol=length(a.vals))
        varpihat <- approx.fn(a.vals,apply(pihat.mat,2,mean),ac)
        varpihat_list[[as.character(c)]] <- varpihat

        muhat.mat <- matrix(as.numeric(stats::predict(mumod,newdata =la.new)$pred),ncol=length(a.vals))
        mahat_list[[as.character(c)]] <- approx.fn(a.vals,apply(muhat.mat,2,mean),ac)
    }
    ylist <- split(y,class_label.factor)
    alist <- split(a,class_label.factor)
    return(drdrtest_em.base(ylist,alist,pihat_list,varpihat_list,
                            muhat_list,mahat_list,arange, h,b,dist, a.grid.size))
}


