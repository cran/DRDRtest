#' The base function for performing tests of average treatment effects
#'
#' This is the base function for testing average treatment effects. Users can use specify the nuisance function values by themselves.
#' 
#' @param y A vector containing the outcomes for each observation
#' @param a A vector containing the treatment levels (dosage) for each observation
#' @param pi A vector containing the propensity scores for each observation
#' @param varpi A vector containing the mean propensity scores for each observation
#' @param mu A vector containing the outcome regression function values for each observation
#' @param ma A vector containing the mean outcome regression fucntion values for each observation
#' @param arange A vector of length 2 giving the lower bound and upper bound of  treatment levels
#' @param h bandwidth to be used in kernel regression. If not specified, will by default use "rule of thumb" bandwidth selector
#' @param b number of Bootstrap samples to be generated
#' @param dist distibution used to generate residuals for Bootstrap samples. Currently only have two options, "TwoPoint" and "Rademachar"
#' @param a.grid.size size of equally spaced grid points over \code{arange} to be generate for numerically evaluating the integral in test statistic
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}
#'     \item{loc.fit:}{A list containg evalution points of average treatment effect and the corresponding values}
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' @export
#' @examples
#' mu.mod<-function(a,l,delta,height){
#'     mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1))+triangle(a-2.5,delta,height)+a*(-0.1*l[,1]+0.1*l[,3])
#'     return(mu)
#' }
#' triangle <- function(a,delta,height){
#'     y <- exp(-a^2/((delta/2)^2))*height
#'     return(y)
#' }
#' set.seed(2000)
#' n <- 500
#' d <- 4
#' sigma <- 0.5
#' delta <- 1
#' height <- 0
#' arange<-c(0.01,4.99)
#' 
#' l <- matrix(rnorm(n*d),ncol=d)
#' colnames(l) <- paste("l",1:4,sep="")
#' logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0.2))
#' lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#' a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#' mu <- mu.mod(a,l,delta,height)
#' residual.list <- rnorm(n,mean=0,sd=sigma)
#' y <- mu+residual.list
#'
#' ## We use the oracal propensity score and outcome regression for illustration
#' pilist <- dbeta(a/5, shape1=lambda, shape2 = 1-lambda)/5
#' varpilist <- colMeans(matrix(dbeta(rep(a,each=n)/5,
#'                                    shape1=rep(lambda,n),
#'                                    shape2 = 1-rep(lambda,n))/5, nrow=n))
#' mulist <- mu
#' malist <-colMeans(matrix(mu.mod(rep(a,each=n),l[rep(1:n,n),],delta,height),nrow=n))
#'
#' out <- drdrtest.base(y,a,pilist,varpilist,mulist,malist,arange)
#' 

drdrtest.base <- function(y,a,pi,varpi,mu,ma,arange,h=NULL,b=1000,dist='TwoPoint',a.grid.size=401){
    pseudo.out <- (y-mu)/(pi/varpi)+ma
    n <- length(y)
    if(is.null(h)){
        h <- 0.9*n^(-1/5)*min(stats::sd(a),stats::IQR(a)/1.34)
    }
    test.out <- wildboot(a,pseudo.out,h,b,arange,dist,a.grid.size)
    return(list(p.value = mean(test.out$T.boot>=test.out$T.obs),
                test.stat = test.out$T.obs,
                Bootstrap.samples = test.out$T.boot,
                loc.fit = test.out$loc.fit,
                bandwidth=h))
}

#' The function for performing tests of average treatment effects with user specified nuisance functions
#'
#' This is the function for testing average treatment effects with user specified nuisance functions.
#' 
#' @param y A vector containing the outcomes for each observation
#' @param a A vector containing the treatment levels (dosage) for each observation
#' @param l A data.frame containing the observations of covariates
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
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}
#'     \item{loc.fit:}{A list containg evalution points of average treatment effect and the corresponding values}
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' @export
#' @examples
#' mu.mod<-function(a,l,delta,height){
#'     mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1))+triangle(a-2.5,delta,height)+a*(-0.1*l[,1]+0.1*l[,3])
#'     return(mu)
#' }
#' triangle <- function(a,delta,height){
#'     y <- exp(-a^2/((delta/2)^2))*height
#'     return(y)
#' }
#' set.seed(2000)
#' n <- 500
#' d <- 4
#' sigma <- 0.05
#' delta <- 1
#' height <- 0
#' arange<-c(0.01,4.99)
#' 
#' l <- matrix(rnorm(n*d),ncol=d)
#' colnames(l) <- paste("l",1:4,sep="")
#' logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0.2))
#' lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#' a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#' mu <- mu.mod(a,l,delta,height)
#' residual.list <- rnorm(n,mean=0,sd=sigma)
#' y <- mu+residual.list
#'
#' ## We use the oracal propensity score and outcome regression for illustration
#' pifunc <- function(a,l){
#'    l <- as.matrix(l)
#'    logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0.2))
#'    lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#'    return(dbeta(a/5,shape1=lambda,shape2 = 1-lambda)/5)
#'}
#'
#' mufunc <- function(a,l){
#'    l <- as.matrix(l)    
#'    return(mu.mod(a,l,delta,height))
#'}
#' out <- drdrtest(y,a,data.frame(l),arange,pifunc,mufunc)

drdrtest <- function(y,a,l,arange,pifunc,mufunc,h=NULL,b=1000,dist='TwoPoint',pi.low=0.01,a.grid.size=401){
    approx.fn <- function(x,y,z){
        stats::predict(stats::smooth.spline(x,y),x=z)$y
    }


    n <- length(y)
    la <- data.frame(l,a=a)
    a.vals <-seq(arange[1],arange[2],length.out=a.grid.size)
    la.new <- rbind(la,data.frame(l[rep(1:n,length(a.vals)),],a=rep(a.vals,each=n)))
    l.new <- la.new[,-ncol(la.new)]
    a.new <- la.new[,ncol(la.new)]

    pihat.vals <- pifunc(a.new, l.new)
    pihat <- pmax(pihat.vals[1:n],pi.low)
    pihat.mat <- matrix(pihat.vals[-(1:n)],nrow = n , ncol= length(a.vals))
    varpihat <- pmax(approx.fn(a.vals,apply(pihat.mat,2,mean),a),pi.low)

    muhat.vals <- mufunc(a.new,l.new)
    muhat <- muhat.vals[1:n]
    muhat.mat <- matrix(muhat.vals[-(1:n)],nrow=n,ncol=length(a.vals))
    mhat <- approx.fn(a.vals,apply(muhat.mat,2,mean),a)

    return(drdrtest.base(y,a,pihat,varpihat,muhat,mhat,arange,h,b,dist,a.grid.size))
    }

#' The function for performing tests of average treatment effects with SuperLearner
#'
#' This is the function for testing average treatment effects with user specified nuisance functions.
#' 
#' @param y A vector containing the outcomes for each observation
#' @param a A vector containing the treatment levels (dosage) for each observation
#' @param l A data.frame containing the observations of covariates
#' @param arange A vector of length 2 giving the lower bound and upper bound of  treatment levels
#' @param pi.sl.lib Models will be used by SuperLearner to estiamte propensity scores
#' @param mu.sl.lib Models will be used by SuperLearner to estiamte outcome regression function
#' @param mu.family Type of response. Currently only support "gaussian" and "binomial"
#' @param h bandwidth to be used in kernel regression. If not specified, will by default use "rule of thumb" bandwidth selector
#' @param b number of Bootstrap samples to be generated
#' @param dist distibution used to generate residuals for Bootstrap samples. Currently only have two options, "TwoPoint" and "Rademachar"
#' @param a.grid.size size of equally spaced grid points over \code{arange} to be generate for numerically evaluating the integral in test statistic
#' @param pi.low Lower bound to truncate propensity scores
#' @param pi.var.low Lower bound to truncate conditional variance of treament (used in propensity score estimation).
#' @return A list containing \describe{
#'     \item{p.value:}{P value of the test result}
#'     \item{test.stat:}{Value of the observed test statistic}
#'     \item{Bootstrap.samples:}{A vector containing test statistic values from Bootstrap samples}
#'     \item{loc.fit:}{A list containg evalution points of average treatment effect and the corresponding values}
#'     \item{bandwidth:}{Bandwidth used in kernel regression}
#' }
#' @import SuperLearner
#' 
#' @export
#' @examples
#' mu.mod<-function(a,l,delta,height){
#'     mu <- as.numeric(l%*%c(0.2,0.2,0.3,-0.1))+triangle(a-2.5,delta,height)+a*(-0.1*l[,1]+0.1*l[,3])
#'     return(mu)
#' }
#' triangle <- function(a,delta,height){
#'     y <- exp(-a^2/((delta/2)^2))*height
#'     return(y)
#' }
#' set.seed(2000)
#' n <- 500
#' d <- 4
#' sigma <- 0.05
#' delta <- 1
#' height <- 0
#' arange<-c(0.01,4.99)
#' 
#' l <- matrix(rnorm(n*d),ncol=d)
#' colnames(l) <- paste("l",1:4,sep="")
#' logit.lambda <- as.numeric(l%*%c(0.1,0.1,-0.1,0.2))
#' lambda <- exp(logit.lambda)/(1+exp(logit.lambda))
#' a <- rbeta(n, shape1 = lambda, shape2 =1-lambda)*5
#'
#' mu <- mu.mod(a,l,delta,height)
#' residual.list <- rnorm(n,mean=0,sd=sigma)
#' y <- mu+residual.list
#'
#' out <- drdrtest.superlearner(y,a,l,arange,pi.sl.lib=c("SL.glm"),mu.sl.lib=c("SL.glm"))

drdrtest.superlearner <- function(y,a,l,arange,pi.sl.lib=c("SL.earth","SL.glm","SL.gam","SL.glmnet"),
                                  mu.sl.lib=c("SL.earth","SL.glm","SL.gam","SL.glmnet"),
                                  mu.family ="gaussian",h=NULL,b=1000,
                                  dist='TwoPoint',a.grid.size=401,pi.low=0.01,pi.var.low=0.01){
    approx.fn <- function(x,y,z){
        stats::predict(stats::smooth.spline(x,y),x=z)$y
    }
    n <- length(y)
    la <- data.frame(l,a=a)
    a.vals <-seq(arange[1],arange[2],length.out=a.grid.size)
    la.new <- rbind(la,data.frame(l[rep(1:n,length(a.vals)),],a=rep(a.vals,each=n)))
    l.new <- la.new[,-ncol(la.new)]

    pimod <- SuperLearner::SuperLearner(Y=a,X=as.data.frame(l),SL.library=pi.sl.lib,newX=l.new)
    pimod.vals <-pimod$SL.predict

    sq.res <- (a-pimod.vals[1:n])^2
    pimod2 <- SuperLearner::SuperLearner(Y=sq.res,X=as.data.frame(l),SL.library=pi.sl.lib,newX=l.new)
    pimod2.vals <- pmax(pimod2$SL.predict,pi.var.low)

    mumod <- SuperLearner::SuperLearner(Y=y,X=data.frame(l,a=a),SL.library=mu.sl.lib,newX=la.new,family=mu.family)
    muhat.vals <- mumod$SL.predict

    a.std <-(la.new$a-pimod.vals)/sqrt(pimod2.vals)
    pihat.vals <- approx.fn(stats::density(a.std[1:n])$x,stats::density(a.std[1:n])$y,a.std)/sqrt(pimod2.vals)
    pihat <- pmax(pihat.vals[1:n],pi.low)

    pihat.mat <- matrix(pihat.vals[-(1:n)],nrow = n , ncol= length(a.vals))
    varpihat <- pmax(approx.fn(a.vals,apply(pihat.mat,2,mean),a),pi.low)

    muhat <- muhat.vals[1:n]
    muhat.mat <- matrix(muhat.vals[-(1:n)],nrow=n,ncol=length(a.vals))
    mhat <- approx.fn(a.vals,apply(muhat.mat,2,mean),a)

    return(drdrtest.base(y,a,pihat,varpihat,muhat,mhat,arange,h,b,dist,a.grid.size))
}
