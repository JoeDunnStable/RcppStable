require("stablecpp")
require("ggplot2")

###--- Tail approximations etc  -- both for pstable() and dstable()
dPareto <- stablecpp:::dPareto

source(system.file("test-tools-1.R", package = "Matrix"), keep.source=interactive())
					#-> identical3(), showProc.time(),...

nc <- 512 # number of points for curve()

pdf("stable-tails.pdf")

pstab.tailratio <- function(alpha, beta, n = nc, prob = 2^-40,
                            xmin = qstable(0.95, alpha,beta, tol = 0.01),
                            xmax = qstable(1 - prob, alpha,beta))
{
  ## Purpose: Explore eps(x) where  (1 - pstable(x))/(1- pPareto()) = 1 + eps(x)
  ##  <==>
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 21 Mar 2011, 10:09
    cl <- match.call()
    stopifnot(0 < prob, prob < .5,
              0 < xmin, xmin < xmax)
    x <- exp(seq(log(xmin), log(xmax), length.out = n))
    iF <-               pstable(x, alpha,beta, lower.tail=FALSE)
    ok <- iF > 0
    iF <- iF[ok]
    x <- x[ok]
    iFp <- stablecpp:::pPareto(x, alpha,beta, lower.tail=FALSE)
    eps <- (iF - iFp)/iFp
    structure(list(x=x, eps=eps, call = cl, alpha=alpha, beta=beta),
              class = "pstableTailratio")
}

plot.pstableTailratio <- function(x, type="l", col="blue3",
                                  lin.col = adjustcolor("red2",.6), ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
    stopifnot(is.list(x), is.numeric(x$eps))
    tit <- substitute("tail ratio approx. for"~~ pstable(alpha==A, beta==B),
                      list(A=x[["alpha"]], B=x[["beta"]]))
    tit<-substitute(atop(num, epsilon(x) == (bar(F)(x,.) - bar(F)[P](x,.)) / bar(F)[P](x,.)),list(num=tit))
    dat <-as.data.frame(x[c("x","eps")])
    dat$type=rep("actual",length(x$x))
#    dat <- dat[dat[,"eps"] > 0, ] ## drop invalid  eps[]
    fm <- lm(log(abs(eps)) ~ log(x), weights = x^2, data = dat)
    dat <- rbind(dat,data.frame(x=dat["x"],eps=exp(predict(fm)),type=rep("fit",nrow(dat))))
    old_theme<-theme_update(legend.position=c(1,1),legend.justification=c(1,1))
    Form <- function(x) formatC(x, digits=4, wid=1)
    fit.line <-substitute(log(epsilon) == A + B * log(x),
                          list(A = Form(coef(fm)[["(Intercept)"]]),
                               B = Form(coef(fm)[["log(x)"]])))
    print(qplot(x=x, y=abs(eps), log = "xy", data = dat, color=type,
                ylab = expression(abs(epsilon) ~~~ "('abs(eps)')"),
                main = tit, geom="line", ...)+
          scale_color_discrete(breaks=c("actual","fit"),
                               labels=c("actual", fit.line))+
          annotation_logticks())
    theme_set(old_theme)
}


plot(tr0  <- pstab.tailratio(1, 0.5))
plot(tr1  <- pstab.tailratio(1.1, 0.25))
plot(tr2 <- pstab.tailratio(0.99, +0.992))

showProc.time()

plot(tr   <- pstab.tailratio(1.2, 0.5))

plot(tr3 <- pstab.tailratio(0.7, +0.9))

plot(tr4 <- pstab.tailratio(1.7, +0.6))

showProc.time()

plot(tr5   <- pstab.tailratio(.1, 0.5))

plot(tr6 <- pstab.tailratio(0.2, +0.9))

plot(tr7 <- pstab.tailratio(.5, +0.6))

showProc.time()

showProc.time()

##---------------- Now the density

##' @title Explore eps(x) where   dstable(x)/dPareto(x) = 1 + eps(x)
##' @param alpha
##' @param beta
##' @param n
##' @param prob
##' @param xmin
##' @param xmax
##' @return an object of \code{\link{class} "dstableTailratio"}, ...
##' @author Martin Maechler, 21 Mar 2011
dstab.tailratio <- function(alpha, beta, n = nc, prob = 2^-40,
                            xmin = qstable(0.95, alpha,beta, tol = 0.01),
                            xmax = qstable(1 - prob, alpha,beta))
{
    cl <- match.call()
    stopifnot(0 < prob, prob < .5,
              0 < xmin, xmin < xmax)
    x <- exp(seq(log(xmin), log(xmax), length.out = n))
    f <-               dstable(x, alpha,beta)
    ok <- f > 0
    f <- f[ok]
    x <- x[ok]
    fp <- stablecpp:::dPareto(x, alpha,beta)
    eps <- (f - fp)/fp
    structure(list(x=x, eps=eps, call = cl, alpha=alpha, beta=beta),
              class = "dstableTailratio")
}

##' @title plot() method for  dstableTailratio() results
##' @param x object of \code{\link{class} "dstableTailratio"}.
##' @param type plot type; default simple lines
##' @param col
##' @param lin.col
##' @param ... optional further arguments passed to \code{\link{plot.formula}()}.
##' @return
##' @author Martin Maechler
plot.dstableTailratio <- function(x, type="l", col="blue3",
                                  lin.col = adjustcolor("red2",.6), ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
    stopifnot(is.list(x), is.numeric(x$eps))
    tit <- substitute("tail ratio approx. for"~~ dstable(alpha==A, beta==B),
                      list(A=x[["alpha"]], B=x[["beta"]]))
    tit<-substitute(atop(num, epsilon(x) == (f(x,.) - f[P](x,.)) / f[P](x,.)),list(num=tit))
    dat <- as.data.frame(x[c("x","eps")])
    dat$type=rep("actual",length(x$x))
    dat <- dat[dat[,"eps"] > 0, ] ## drop invalid  eps[]
    fm <- lm(log(eps) ~ log(x), weights = x^2, data = dat)
    dat <- rbind(dat,data.frame(x=dat["x"],eps=exp(predict(fm)),type=rep("fit",nrow(dat))))
    ## NOTA BENE: Empirically,  I see that  eps > 0 <==> alpha > 1
    ##                                      eps < 0 <==> alpha < 1
    Form <- function(x) formatC(x, digits=4, wid=1)
    fit.line <-substitute(log(epsilon) == A + B * log(x),
                          list(A = Form(coef(fm)[["(Intercept)"]]),
                               B = Form(coef(fm)[["log(x)"]])))
    gph<-qplot(x=x, y=eps, log = "xy", data = dat, color=type,
                ylab = expression(epsilon ~~~ "('eps')"),
                main = tit, geom="line", ...)+
            scale_color_discrete(breaks=c("actual","fit"),
                                 labels=c("actual", fit.line))+
            annotation_logticks()
    gph
}

old_theme<-theme_update(legend.position=c(1,1),legend.justification=c(1,1))

plot(fr   <- dstab.tailratio(1.01, 0.8))
plot(fr   <- dstab.tailratio(1.05, 0.4))
plot(fr   <- dstab.tailratio(1.1,  0.4))
plot(fr   <- dstab.tailratio(1.2,  0.5))
plot(fr   <- dstab.tailratio(1.3,  0.6))

showProc.time()

plot(fr   <- dstab.tailratio(1.4, 0.7))
plot(fr   <- dstab.tailratio(1.5, 0.8))

plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1000))
plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1e4))+geom_vline(aes(xintercept=1000), linetype=2)
plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1e5))+geom_vline(aes(xintercept=1e4), linetype=2)

showProc.time()

plot(fr   <- dstab.tailratio(1.6, 0.9))
plot(fr   <- dstab.tailratio(1.7, 0.1))
plot(fr   <- dstab.tailratio(1.8, 0.2))
theme_set(old_theme)

showProc.time()

##------ Some explicit tail problems visualized:

I <- integrate(dstable, 0,Inf, alpha=0.998, beta=0, subdivisions=1000)
str(I)
stopifnot(abs(I$value - 0.5) < 1e-4)

qcurve<-function(alpha, beta, xmin, xmax, pareto_flag=F, log_flag=T){
  x<-exp(seq(log(xmin),log(xmax),length.out=nc))
  y<-dstable(x, alpha=alpha, beta=beta,    log=T)
  df<-data.frame(x=x,y=y,type=rep("dstable", nc))
  if (pareto_flag){
    y2<-dPareto(x, alpha=alpha, beta=beta, log=T)
    df<-rbind(df,data.frame(x=x,y=y2,type=rep("dPareto",nc)))
  }
  q_log<-if (log_flag) "x" else ""
  s_log<-if (log_flag) "b" else ""
  print(qplot(x=x,y=y, color=type, data=df, geom="line", log=q_log,
      main=bquote(dstable(x, alpha == .(alpha), beta == .(beta), log == T)))+
        annotation_logticks(sides=s_log))
}

qcurve(alpha=.999, beta=0.1, xmin=10, xmax=1e17)
qcurve(alpha=.999, beta=0.9, xmin=10, xmax=1e17)
qcurve(alpha=.999, beta=0.99, xmin=10, xmax=1e17)
qcurve(alpha=.999, beta=0.99, xmin=10, xmax=1e170)
showProc.time()

## less problems when alpha > ~= 1
qcurve(alpha=1.001, beta=0.99, xmin=10, xmax=1e7)
qcurve(alpha=1.001, beta=0.99, xmin=10, xmax=1e17)
## -> zoom in on former problem:
qcurve(alpha=1.001, beta=0.99, xmin=1e12, xmax=160e12, log=F, pareto_flag=T)

qcurve(alpha=1.001,beta=0.99, xmin=10, xmax=1e40)
showProc.time()

## NB: alpha == 1   also has been fixed
qcurve(alpha=1., beta=0.99, xmin=1, xmax=20, log_flag=F)
qcurve(alpha=1., beta=0.99, xmin=1, xmax=100, log_flag=F)

dev.off()
showProc.time()
