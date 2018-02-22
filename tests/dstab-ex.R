require("RcppStable")
require("ggplot2")
require("plyr")

dPareto <- RcppStable:::dPareto

source(system.file("test-tools-1.R", package = "Matrix"), keep.source=interactive())
					#-> identical3(), showProc.time(),...

stopifnot(0 < print(dstable(4000., alpha=1.00001, beta=0.6)))
## gave error in fBasics::dstable()
## now 18 warnings from uniroot()

pdf("dstab-ex.pdf")

x <- 2^seq(0, 20, length= 200)
## Regression check for alpha=2: <==> *norm() :
x. <- x/1024
fx <- dstable(x., alpha = 2, beta = 0, gamma = 1.1, delta=2, pm=2)
lf <- dstable(x., alpha = 2, beta = 0, gamma = 1.1, delta=2, pm=2, log=TRUE)
stopifnot(
    local({i <- is.finite(log(fx)); all.equal(log(fx[i]), lf[i])}),
    all.equal(fx, dnorm(x., m=2, s=1.1)),
    all.equal(lf, dnorm(x., m=2, s=1.1, log=TRUE)))

fx <- dstable(x, alpha = 1.0001, beta = 0.6)

qplot(x,fx, log="x", geom="line",
      main = expression(dstable(x, alpha == 1.0001, beta == 0.6)))# looks good
qplot(x,fx, log="xy", geom="line",
      main = expression(dstable(x, alpha == 1.0001, beta == 0.6)))# --- perfect now
stopifnot((dlx <- diff(log(fx))) < 0, # decreasing
      { dl <- dlx[x[-1] > 1000] # the negative slope:
	relErr(-0.13934, dl) < 4e-4},
	  diff(dl) > -1e-6)#  > 0 "in theory"; > -6e-7, ok on 64-bit Linux

nc <- 512 # number of points for curve()

zeta <- function(alpha,beta) if(alpha==1) 0 else -beta*tan(pi/2*alpha)

## negative beta:
x=seq(-.5,1.5,length.out=nc)
qplot(x=x,y=dstable(x,.75,-.5),ylab="f(x)",main=expression(dstable(x, alpha == .75, beta == -.5)),geom="line")
## cx <- curve(dstable(x, 0.75, -.5), -.5, 1.5, n=nc)# ok, now
m <- stableMode(0.75, -.5, tol=1e-14)
stopifnot(all.equal(m, 0.35810298366, tol = 1e-7),
	  all.equal(dstable(m, 0.75, -.5), 0.3554664043, tol=1e-6))

showProc.time()

###-------- "small" alpha -----------------
## alpha --> 0 --- very heavy tailed -- and numerically challenging.

## symmetric (beta = 0)
(x0 <- (-16:16)/256)
fx0 <- dstable(x0, alpha = 0.1, beta=0, gamma = 1e6)
qplot(x=x0, y=fx0, geom=c("line"),
     main = expression(f(x, alpha== 0.1, beta == 0, gamma == 10^6)))+
  geom_point(shape="o",size=5)
stopifnot(all.equal(fx0[17],1.15508291498374),
	  all.equal(fx0[ 1],0.02910420736536),
	  all.equal(range(diff(fx0[1:8])),
		    c(0.0011871409, 0.0025179435), tol=1e-6)
	  )

## beta > 0
x<-seq(-1,1,length.out=nc)
qplot(x=x,y=dstable(x, alpha = 0.3, beta = 0.5, tol=1e-7),ylab="f(x)",
            main=expression(dstable(x, alpha== 0.3, beta  == 0.5, tol == 10^-7)),geom="line")
m3 <- stableMode(0.3, 0.5, tol=1e-14)# still with 3 warnings
stopifnot(all.equal(m3, -0.2505743952946, tol = 2e-8))
## zoom into the above
x<-seq(-.27, -.22, length.out=nc)
qplot(x=x,y=dstable(x, alpha = 0.3, beta = 0.5, tol=1e-7), ylab="f(x)",
      main=expression(dstable(x, alpha== 0.3, beta  == 0.5, tol == 10^-7)),geom="line")+
    geom_vline(aes(xintercept = m3), color="green", linetype=2, size=1)

x<-seq(-.4,.2,length.out=nc)
y<-dstable(x, alpha = 0.1, beta = 0.5, tol=1e-7)
m1<-stableMode(.1,.5,tol=1e-15)
qplot(x=x,y=y,geom="line",
      main=expression(dstable(x, alpha == 0.1, beta == 0.5, tol == 10^-7)),
      ylab = "f(x)", ylim = c(0, 25))+
  geom_vline(aes(xintercept=m1),color="green",linetype=2,size=1)+
  geom_hline(aes(yintercept=0), color="red", linetype=2,size=1)
stopifnot(all.equal(m1, -0.079192219, tol=2e-7)) # -0.079192219, was -0.0791921758
## check mode *and* unimodality
i. <- x > m1
stopifnot(## decreasing to the right:
	  diff(y[ i.]) < 0,
	  ## increasing on the left:
	  diff(y[!i.]) > 0)

## Now *really* small alpha --> 0:
##     --------------------------
## Note that
if(require("FMStable")) {
    try( FMStable::setParam(alpha = 1e-18, location=0, logscale=0, pm = 0) )
## gives
## Error in .... setParam: S0=M parametrization not suitable for small alpha
}
## now if this is true (and I think we can trust Geoff Robinson on this),
## we currently have a  "design bug - problem":
## as internally, we always translate to  'pm=0' parametrization __FIXME_(how ??)__
## --> solution: see below:  there's a simple  (alpha -> 0) asymptotic

## These now "work": .... well with integration warnings !!

pdstab.alpha <- function(x, betas, alphas = 2^-(40:2),title, ...)
{
    stopifnot(is.numeric(x),
              is.numeric(betas), -1 <= betas, betas <= 1,
              is.numeric(alphas), 0 <= alphas, alphas <= 2)
    if(is.unsorted(alphas)) alphas <- sort(alphas)
    if(is.unsorted(beta)) betas<-sort(betas)
    one_ab<-function(df){
      data.frame(x=x,density=dstable(x, alpha=df$alpha, beta = df$beta, pm = 0))
    }
    df<-expand.grid(alpha=alphas,beta=betas)
    df_out<-ddply(df,.(alpha,beta),one_ab)
    qplot(x=alpha,y=density,data=df_out,color=as.factor(beta),
             linetype=as.factor(x), geom="line",log="xy",
             xlab = quote(alpha),
             ylab = expression(f(x, alpha)),
             main = title)+
      scale_color_discrete(guide_legend(title="beta"))+
      scale_linetype_discrete(guide_legend(title="x"))
}

## vary beta,  keep x :
pdstab.alpha(x = -1, betas = c(0.1,.3,.5,.7,.9),
             title=expression(dstable(x == -1, alpha, beta, pm == 0)))
## vary x,  keep beta
pdstab.alpha(x =  c(-1,+1,+5,+50,-10), beta = 0.3,
             title=expression(dstable(x, beta == .3, alpha, pm == 0)))

## Plots suggest a simple formula
##  log(f(x, alpha, beta))=  c(x,beta) + log(alpha)
##      f(x, alpha, beta) =  C(x,beta) * (alpha)   -- ???

## for normal alpha, it looks different {which is reassuring!}
pdstab.alpha(x = -1, betas = 0.3, alphas = seq(1/128, 2, length=100),
             title=expression(dstable(x == -1, alpha, beta, pm == 0)))

## Show the formula, we derived:
## f(x, \alpha, \beta) ~  \alpha * \frac{1 + \beta}{2e \abs{x + \pi/2 \alpha\beta}}
## ONLY ok, when  x > zeta := - beta * tan(pi/2 *alpha)
## otherwise, "swap sign" of (x, beta)
dst.smlA <- function(x, alpha, beta, log = FALSE) {
    pa <- pi/2*alpha
    i. <- x < -pa*beta
    if(any(i.)) {
        beta <- rep(beta, length.out=length(x))
        beta[i.] <- -beta[i.]
        x   [i.] <- -x   [i.]
    }
    ## alpha*(1+beta) / (2*exp(1)*(x+ pa*beta))
    r <- log(alpha) + log1p(beta) - (1 + log(2*(x+ pa*beta)))
    if(log) r else exp(r)
}

set.seed(17)

alpha <- 1e-15 ## must be larger than cutoff in .fct1() in ../R/dpqr-stable.R :
for(beta in runif(20, -1, 1)) {
 cat(sprintf("beta =%10g: ", beta))
 for(N in 1:10) {
     x <- 10*rnorm(100); cat(".")
     stopifnot(all.equal(dstable (x, alpha, beta),
                         dst.smlA(x, alpha, beta), tol = 1e-13))
 }; cat("\n")
}


x=seq(-10,10,length.out=101)
df<-rbind(data.frame(x=x,y=dstable (x, alpha=1e-10, beta=.8, log=TRUE),type=rep("dstable",101)),
          data.frame(x=x,y=dst.smlA(x, alpha=1e-10, beta=.8, log=TRUE),type=rep("dst.smlA",101)))
qplot(x=x,y=y,data=df,color=type,geom="line",
          main=expression(f(x, alpha == 10^-10, beta == .8)))

## "same" picture (self-similar !)
x<-seq(-1,1,length.out=101)
df<-rbind(data.frame(x=x,y=dstable (x, alpha=1e-10, beta=.8, log=TRUE),type=rep("dstable",101)),
          data.frame(x=x,y=dst.smlA(x, alpha=1e-10, beta=.8, log=TRUE),type=rep("dst.smlA",101)))
qplot(x=x,y=y,data=df,color=type,geom="line",
      main=expression(f(x, alpha == 10^-10, beta == .8)))



## Testing stableMode here:


### beta = 1 (extremely skewed)  and small alpha: ---------
##  --------
## Problem at *left* ("less problematic") tail, namely very close to where the
## support of the density becomes mathematically exactly zero :
##
## clear inaccuracy / bug --- maybe *seems* curable
##
x<-seq(-40,40,length.out=101)
qplot(x=x,y=dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE),geom="line",
      ylab="y",main=expression(dstable(e^x, alpha == 0.1, beta == 1, pm == 1, log == TRUE)))
##            ------
## --> warnings both from uniroot ("-Inf") and .integrate2()
## about equivalent to
x<-exp(seq(log(1e-15),log(4e4),length.out=1001))
qplot(x=x,y=dstable(x, alpha= 0.1, beta=1, pm=1, log=TRUE), geom="line",log="x",
      main=expression(dstable(e^x, alpha == 0.1, beta == 1, pm == 1, log == TRUE)))

x<-seq(-40,20,length.out=101)
qplot(x=x,y=dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE), geom="line",
      main=expression(dstable(e^x, alpha == 0.1, beta == 1, pm == 1, log == T)))
## or here, ... but still not good enough
x<-seq(-45,30,length.out=101)
qplot(x=x,y=dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE), geom="line",
      main=expression(dstable(e^x, alpha == 0.1, beta == 1, pm == 1, log == T)))

showProc.time()

##------ NB: Pareto tail behavior -- see more in ./tails.R
##						   =======

## alpha ~= 1  ---- and x ~ zeta(a,b) -----------
## ==========
f1 <- dstable(6366.197,	 alpha= 1.00001, beta= .1)
f2 <- dstable(-50929.58, alpha= 1.00001, beta= -.8)
stopifnot(f1 > 0, f2 > 0)

## these all work (luck):
chk_near_zeta<-function(alpha,beta){
  zet <- zeta(alpha= alpha, beta= beta)
  d0<-dstable(zet,alpha,beta)
  ## here, we must have larger zeta.tol: = 5e-5 is fine; now using "automatic" default
  x<-max(abs(zet),1e-8/abs(1-alpha))*seq(-5e-4,5e-4,length.out=nc+if(nc%%2==0) 1 else 0)
  df<-data.frame(x=x,y=dstable(zet+x, alpha= alpha, beta= beta),fnct="dstable")
  if (abs(zet)>1e4)
    df<-rbind(df,data.frame(x=x,y=dPareto(zet+x, alpha= alpha, beta= beta),fnct="dPareto"))
  z.txt <- bquote(paste(x == zeta(.), phantom() == .(signif(zet,6))))
  show(qplot(x=x,y=y, geom="line",data=df,color=fnct,ylim=c(0,max(df$y)),
        main=bquote(dstable(x+zeta(alpha,beta), alpha == .(alpha), beta == .(beta)))) +
        geom_hline(aes(yintercept=d0),color="green")+
        geom_vline(aes(xintercept=0), color="pink")+
        geom_text(aes(x=0,y=max(df$y)-.1*(max(df$y)-min(df$y)),label=deparse(z.txt)), alpha=1,parse=T,hjust=0, color="pink"))
  df
}
## no longer much noise (thanks to zeta.tol = 1e-5):
df<-chk_near_zeta(1.00001,-.8)
stopifnot({ rr <- range(df$y)
	    2.1e-10 < rr & rr < 2.3e-10 })
df<-chk_near_zeta(1.00001,0)

df<-chk_near_zeta(.1,0)

showProc.time()

### ---- alpha == 1 ---------
x<-seq(-20,20,length.out=nc)
qplot(x=x,y=dstable(x, alpha = 1, beta = 0.3), geom="line", log="y",
      main=bquote(dstable(x, alpha == 1, beta == 0.3)),ylab="f(x)")

x<-seq(-200, 160, length.out=nc)
df<-rbind(data.frame(x=x,y=dstable(x, alpha = 1, beta = 0.3, log=TRUE),fnct=rep("dstable",nc)),
          data.frame(x=x,y=dPareto(x, alpha = 1, beta = 0.3, log=TRUE),fnct=rep("dPareto",nc)))
qplot(x=x,y=y,data=df,color=fnct, geom="line",
      main=bquote(dstable(x, alpha == 1, beta == 0.3, log == T)),
      ylab="f(x)")
## was discontinuous but no longer
## ditto:
x<-seq(-70,80,length.out=nc)
df<-rbind(data.frame(x=x,y=dstable(x, alpha=1, beta= 0.1, log=TRUE),fnct="dstable"),
          data.frame(x=x,y=dPareto(x, alpha=1, beta= 0.1, log=TRUE),fnct="dPareto"))
qplot(x=x, y=y, data=df, color=fnct, geom="line",
      main=bquote(dstable(x, alpha == 1, beta == 0.1, log == T)),
      ylab="f(x)")

showProc.time()

dstable(-44, alpha=1, beta= .1)# failed
## large x gave problems at times:
dstable(-1e20, alpha = 0.9,  beta = 0.8)

chkUnimodal <- function(x) {
    ## x = c(x1, x2)  and  x1 is *increasing*  and x2 is *decreasing*
    stopifnot((n <- length(x)) %% 2 == 0,
	      (n2 <- n %/% 2) >= 2)
    if(is.unsorted(x[seq_len(n2)])) warning("first part is *not* increasing", immediate.=T)
    if(is.unsorted(x[n:(n2+1)]))    warning("seconds part is *not* decreasing", immediate.=T)
    invisible(x)
}

showProc.time()

xLrg <- c(10^c(10:100,120, 150, 200, 300), Inf)
xLrg <- sort(c(-xLrg, xLrg))
d <- dstable(xLrg, alpha = 1.8,	  beta = 0.3 ); chkUnimodal(d)
d <- dstable(xLrg, alpha = 1.01,  beta = 0.3 ); chkUnimodal(d) # (was slow!)
## look at the problem:
dstCurve <- function(alpha, beta, log=TRUE, NEG=FALSE,
		     from, to, n=nc, cLog="", ...)
{
  x<-exp(seq(from=log(from),to=log(to),length.out=n))
  log_flag<-if (log) "TRUE" else "FALSE"
  if(NEG) {
    df<<-rbind(data.frame(x=x,y=dstable(-x, alpha=alpha, beta=beta, log=log),fnct=rep("dstable",nc)),
             data.frame(x=x,y=dPareto(-x, alpha=alpha, beta=beta, log=log),fnct=rep("dPareto",nc)))
    qplot(x=x,y=y, data=df,geom="line", color=fnct, log=cLog,
          main=bquote(dstable(-x, alpha == .(alpha), beta == .(beta), log == .(log_flag))),
          ylab="f(x)", ...)
  } else {
    df<<-rbind(data.frame(x=x,y=dstable(x, alpha=alpha, beta=beta, log=log),fnct=rep("dstable",nc)),
             data.frame(x=x,y=dPareto(x, alpha=alpha, beta=beta, log=log),fnct=rep("dPareto",nc)))
    qplot(x=x,y=y, data=df, geom="line", color=fnct, log=cLog,
          main=bquote(dstable(x, alpha == .(alpha), beta == .(beta), log == .(log_flag))),
          ylab="f(x)",...)
  }
}

## dstable and dPareto are slightly different at this scale
dstCurve(alpha = 1.01, beta = 0.3, NEG=TRUE, log=FALSE, from=1e10, to=1e40 ,cLog="xy")
dstCurve(alpha = 1.01, beta = 0.3, NEG=FALSE, log=FALSE, from=1e10, to=1e40 ,cLog="xy")
## zoom in:
dstCurve(alpha = 1.01, beta = 0.3, NEG=TRUE, log=FALSE, from=1e15,to= 1e16,cLog="xy")
dstCurve(alpha = 1.01, beta = 0.3, NEG=TRUE, log=FALSE, from=1e39,to= 1e40,cLog="xy")
dstCurve(alpha = 1.01, beta = 0.3, NEG=FALSE, log=FALSE, from=1e15,to= 1e16,cLog="xy")
dstCurve(alpha = 1.01, beta = 0.3, NEG=FALSE, log=FALSE, from=1e39,to= 1e40,cLog="xy")
showProc.time()

d <- dstable(xLrg, alpha = 1.001, beta = -0.9) # No warnings (had over 50)
chkUnimodal(d)
## look at the former localtion of the problem
dstCurve(alpha = 1.001, beta = -0.9, log=FALSE, NEG=TRUE,  1e10, 1e40, cLog="xy")
## and at the right tail, too:
dstCurve(alpha = 1.001, beta = -0.9, log=FALSE, NEG=FALSE, 1e10, 1e40, cLog="xy")

d <- dstable(xLrg, alpha = 1. ,	  beta = 0.3 ); chkUnimodal(d) # "ok" now
d <- dstable(xLrg, alpha = 0.9,	  beta = 0.3 ) # No warnings (had 11)
chkUnimodal(d)
d <- dstable(xLrg, alpha = 0.5,	  beta = 0.3 ) # No warnings (had 22)
chkUnimodal(d)
d <- dstable(xLrg, alpha = 0.1,	  beta = 0.3 ) # Over 50 warnings (had 21)
chkUnimodal(d)
## look at the problem:
dstCurve(alpha = .1, beta = .3, log=FALSE, NEG=TRUE,  1e10, 1e40, cLog="xy")
## Zoom in
dstCurve(alpha = .1, beta = .3, log=FALSE, NEG=TRUE,  1e39, 1e40, cLog="xy")
## and at the right tail, too:
dstCurve(alpha = .1, beta = .3, log=FALSE, NEG=FALSE, 1e10, 1e40, cLog="xy")
## Zoom in
dstCurve(alpha = .1, beta = .3, log=FALSE, NEG=FALSE,  1e39, 1e40, cLog="xy")

showProc.time()

##-------------	 beta = 1  ---------------------
options(dstable.debug = TRUE)
dstable(1, alpha=1.2,	beta= 1 - 1e-7)#ok
dstable(1, alpha=1.2,	beta= 1)# gave error, because	g(pi/2) < 0
dstable(0, alpha=13/16, beta= 1 -2^-52)# had error as	g(-theta0)  |->	 NaN
dstable(0, alpha=19/16, beta= 1)       # had error as	g(pi/2)	    |->	 NaN
options(dstable.debug = FALSE)

## NB: (beta=1, alpha = 1/2) is 'Levy' ---> dLevy() and some checks
## -- in ./pstab-ex.R
##	   ~~~~~~~~~~

require(plyr)
ep <- 2^-(1:54)## beta := 1 - ep ---> 1  {NB: 1 - 2^-54 == 1  numerically}
alph.s <- (1:32)/16   # in (0, 2]
df<-expand.grid(alpha=alph.s,ep=ep)
df_b1<-ddply(df,.(alpha,ep),function(df){data.frame(x=0:10,y=dstable(0:10,df$alpha,1-df$ep))})
print(summary(df_b1$y))
r.b1 <- range(df_b1$y)
stopifnot(0 < r.b1[1], r.b1[2] < 0.35)
## "FIXME" test more: monotonicity in x {mode is < 0 = min{x_i}}, beta, alpha, ...
showProc.time()

beta<-(-4:4)/4
x <- (-80:80)/8

graph_pm2<-function(alpha, betas, x) {
  df_in<-expand.grid(alpha=alpha, beta=beta)

  show_pm2 <- function(df_in) {
    data.frame(x=x, pdf=dstable(x, df_in$alpha, df_in$beta, gamma=3, delta=-1,pm=2))
  }

  df_out <- ddply(df_in, .(alpha, beta), show_pm2)
  df_out$beta <- as.factor(df_out$beta)
  df_vline = data.frame(xintercept=-1)
  tmp <- paste("dstable(x, ", expression(  alpha ), " = ", eval(alpha), ", beta, gamma = 3, delta = -1, pm = 2)")
  gph <- qplot(x=x, y=pdf, data=df_out, color=beta, geom="line",
               main=tmp)+
    geom_vline(aes(xintercept=xintercept),data=df_vline)
  print(gph)
}

graph_pm2(.1, beta, x)
graph_pm2(.5, beta, x)
graph_pm2(.9999, beta, x)
graph_pm2(1, beta, x)
graph_pm2(1.0001, beta, x)
graph_pm2(1.5, beta, x)
graph_pm2(1.9, beta, x)

showProc.time()

alpha<-c(.5, 1, 1.5)
beta<-.5
gamma<-3
delta<--5
pm<-c(0,1,2)

df_r<-expand.grid(alpha=alpha, beta=beta, gamma=gamma, delta=-5, pm=pm)

ks_test <- function(df) {
  set.seed(200)
  tmp<-ks.test(rstable(10000, alpha=df$alpha, beta=df$beta,
                             gamma=df$gamma, delta=df$delta, pm=df$pm),
               "pstable", alpha=df$alpha, beta=df$beta,
                          gamma=df$gamma, delta=df$delta, pm=df$pm)
  data.frame(D=tmp$statistic, p.value=tmp$p.value)
}

df_r_out <- ddply(df_r, .(alpha, beta, gamma, delta, pm), ks_test)
df_r_out

stopifnot(min(df_r_out$p.value)>.01)
showProc.time()

dev.off()
