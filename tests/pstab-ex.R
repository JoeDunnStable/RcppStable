require("stablecpp")
require("ggplot2")
require("grid")

pPareto <- stablecpp:::pPareto

source(system.file("test-tools-1.R", package = "Matrix"), keep.source=interactive())
					#-> identical3(), showProc.time(),...

options(pstable.debug = FALSE)
# options(pstable.debug = TRUE)# want to see when uniroot() is called

stopifnot(all.equal(pstable(0.3, 0.75, -.5),
		    0.66887227658457, tol = 1e-10))
## was                    0.6688726496, tol = 1e-8))

pstable(-4.5, alpha = 1, beta = 0.01)## gave integration error (now uniroot..)

## a "outer vectorized" version:
pstabALL <- function(x, alpha, beta, ...)
    sapply(alpha, function(alph)
	   sapply(beta, function(bet)
			pstable(x, alph, bet, ...)))

alph.s <- (1:32)/16   # in (0, 2]
beta.s <- (-16:16)/16 # in [-1, 1]
stopifnot(pstabALL( Inf, alph.s, beta.s) == 1,
	  pstabALL(-Inf, alph.s, beta.s, log.p=TRUE) == -Inf,
	  pstabALL( 0,	 alph.s, beta = 0) == 0.5,
	  TRUE)

pdf("pstab-ex.pdf")

##---- log-scale -------------
x<-seq(5, 150, length.out=500)
df<-rbind(data.frame(x=x,y=pstable(x, alpha=1.8, beta=.9, lower.tail=FALSE, log.p=TRUE),fnct="pstable"),
          data.frame(x=x,y=pPareto(x, alpha=1.8, beta=.9, lower.tail=FALSE, log.p=TRUE),fnct="pPareto"))
qplot(x=x,y=y,data=df,color=fnct, ylab="log(p(x))",geom="line", log="x",
      main=bquote(pstable(x, alpha == 1.8, beta == .9, lower.tail == F, log.p == T)))+
      annotation_logticks(sides="b")
##--> fixed from R version via use of regions

## the less extreme part - of that:
x<-seq(1, 50, length.out=500)
df<-rbind(data.frame(x=x,y=pstable(x, alpha=1.8, beta=.9, lower.tail=FALSE, log.p=TRUE),fnct="pstable"),
          data.frame(x=x,y=pPareto(x, alpha=1.8, beta=.9, lower.tail=FALSE, log.p=TRUE),fnct="pPareto"))
qplot(x=x,y=y,data=df,color=fnct, ylab="log(f(x))",geom="line", log="x",
      main=bquote(pstable(x, alpha == 1.8, beta == .9, lower.tail == F, log.p == T)))+
      annotation_logticks(sides="b")

## Check that	pstable() is the integral of dstable() --- using simple Simpson's rule

vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

chk.pd.stable <- function(alpha, beta, xmin=NA, xmax=NA,
			  n = 256, do.plot=TRUE,
			  comp.tol = 1e-13, eq.tol = 1e-10)
{
    stopifnot(n >= 20)
    if(is.na(xmin)) xmin <- qstable(0.01, alpha, beta)
    if(is.na(xmax)) xmax <- qstable(0.99, alpha, beta)
    dx <- ceiling(1024*grDevices::extendrange(r = c(xmin, xmax), f = 0.01))/1024
    h <- diff(dx)/n
    x <- seq(dx[1], dx[2], by = h)
    int_dstable<-function(i) {
      integrate(dstable,lower=x[i-1],upper=x[i],alpha=alpha,beta=beta,tol=comp.tol)$value
    }
    Fx <- pstable(x, alpha=alpha, beta=beta, tol=2*comp.tol)
    ## integrate from x[1] up to x[i]
    ## the exact value will be F(x[i]) - F(x[1]) == Fx[i] - Fx[1]
    Fx.<- cumsum(c(0,vapply(seq_along(x)[-1],int_dstable,0.)))
    a.eq <- all.equal(Fx., Fx - Fx[1], tol = eq.tol)
    if(do.plot) {
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(2,1)))
    	## Show the fit
      df<-rbind(data.frame(x=x,y=Fx-Fx[1],fnct="pstable"),
                data.frame(x=x,y=Fx.,fnct="int(dstable"))
      old<-theme_update(legend.position=c(1,0),legend.justification=c(1,0))
    	print(qplot(x=x, y=y, color=fnct,data=df,geom="line",
    	      main=bquote(pstable(x, alpha == .(alpha), beta == .(beta))),
    	      ylab=sprintf("Integral from %g to x",df$x[1])),vp=vplayout(1,1))
    	theme_set(old)
    	## show the "residual", i.e., the relative error
      print(qplot(x=x, y=1- Fx./(Fx - Fx[1]),
    	     geom="line", xlim = range(x),
    	     main="Relative error of pstable vs integral of dstable"),vp=vplayout(2,1))
    }

    if(!isTRUE(a.eq)) stop(a.eq)
    invisible(list(x=x, F=Fx, i. = seq_along(x), F.appr. = Fx.))
}

op <- par(mfrow=2:1, mar = .1+c(3,3,1,1), mgp=c(1.5, 0.6,0))

c1 <- chk.pd.stable(.75, -.5,  -1, 1.5, eq.tol = .006)
c2 <- chk.pd.stable(.95, +0.6, -1, 1.5, eq.tol = .006)# with >= 50 warnings
## here are the "values"
with(c1, all.equal(F.appr., F[i.] - F[1], tol = 0)) # (.0041290 on 64-Lnx)
with(c2, all.equal(F.appr., F[i.] - F[1], tol = 0)) # (.0049307 on 64-Lnx)

showProc.time() #

c3 <- chk.pd.stable(.95, +0.9, -3, 15) # >= 50 warnings

grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
x<-seq(-20, 5,length.out=101)
print(qplot(x=x,y=dstable(x, .999, -0.9),  log="y",geom="line",
      main=bquote(dstable(x, alpha == .999, beta == -.9))),vp=vplayout(1,1))
print(qplot(x=x,y=pstable(x, .999, -0.9), log="y", geom="line",
      main=bquote(pstable(x, alpha == .999, beta == -.9))),vp=vplayout(2,1))#-> using uniroot
c4 <- chk.pd.stable(.999, -0.9, -20, 5)

showProc.time() #

## alpha == 1 , small beta  ---- now perfect
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
x<-seq(-6,9,length.out=101)
print(qplot(x=x, y=pstable(x, alpha=1, beta= .01), geom="line",ylim=0:1,
      main=bquote(pstable(x, alpha == 1, beta == .01)))+
      geom_hline(aes(yintercept=0:1),linetype=3, color="grey")+
      geom_vline(aes(xintercept=0),linetype=3, color="grey"),vp=vplayout(1,1))
n <- length(x <- seq(-6,8, by = 1/16))
px <- pstable(x, alpha=1, beta= .01)
## now take approximation derivative by difference ratio:
## and compare with dstable() ... which used to show dstable() problem now fixed:
x. <- (x[-n]+x[-1])/2
fx. <- dstable(x., alpha=1, beta=.01)
df<-rbind(data.frame(x=x.,y=diff(px)*16,fnct=rep("diff(px)",n-1)),
          data.frame(x=x.,y=fx.,fnct=rep("dstable",n-1)))
old_theme<-theme_update(legend.position=c(1,1),legend.justification=c(1,1))
print(qplot (x=x, y=y, data=df, color=fnct, geom="line",
       main="Approximate derivative of px vs dstable"),vp=vplayout(2,1))
theme_set(old_theme)
## now check convexity/concavity :
f2 <- diff(diff(px))
stopifnot(f2[x[-c(1,n)] < 0] > 0,
	  f2[x[-c(1,n)] > 0] < 0)
x<-seq(-6,50,length.out=500)
qplot(x=x,y=dstable(x, 1., 0.99), geom="line", log="y",
      main=bquote(dstable(x, alpha == 1., beta == .99)), ylab="y")# used to br "uneven" (x < 0); 50 warnings
x<-seq(-10, 30, length.out=500)
qplot(x=x,y=dstable(x, 1.001, 0.95), geom="line",log="y",
      main=bquote(dstable(x, alpha == 1.001, beta == .95)), ylab="y")# much better

showProc.time() #

c5 <- chk.pd.stable(1.,	   0.99,  -6, 50)
c6 <- chk.pd.stable(1.001, 0.95, -10, 30)
with(c5, all.equal(F.appr., F[i.] - F[1], tol = 0)) # .00058 on 64-Lnx
with(c6, all.equal(F.appr., F[i.] - F[1], tol = 0)) # 2.611e-5 on 64-Lnx

## right tail:
c1.0 <- chk.pd.stable(1., 0.8,	-6, 500,n=1000)# uniroot; rel.difference = .030

## show it more clearly
x<-seq(20, 800, length.out=500)
df<-rbind(data.frame(x=x,y=pstable(x, alpha=1, beta=0.5),fnct=rep("pstable",500)),
          data.frame(x=x,y=pPareto(x, alpha=1, beta=0.5),fnct=rep("pPareto",500)))
qplot(x=x, y=y, color=fnct, data=df, geom="line",  log="x", ylim=c(.97, 1),
      main=bquote(pstable(x, alpha == 1, beta == 0.5)))+
  annotation_logticks(sides="b")
# and similarly (alpha ~= 1, instead of == 1):
x<-seq(20, 800,length.out=500)
df<-rbind(data.frame(x=x,y=pstable(x, alpha=1.001, beta=0.5),fnct=rep("pstable",500)),
          data.frame(x=x,y=pPareto(x, alpha=1.001, beta=0.5),fnct=rep("pPareto",500)))
qplot(x=x, y=y, color=fnct, data=df, geom="line",  log="x", ylim=c(.97, 1),
      main=bquote(pstable(x, alpha == 1.001, beta == 0.5)))+
  annotation_logticks(sides="b")
## zoom
x<-seq(100, 200,length.out=500)
df<-rbind(data.frame(x=x,y=pstable(x, alpha=1.001, beta=0.5),fnct=rep("pstable",500)),
          data.frame(x=x,y=pPareto(x, alpha=1.001, beta=0.5),fnct=rep("pPareto",500)),
          data.frame(x=x,y=pstable(x, alpha=1, beta=0.5),fnct=rep(deparse(bquote(pstable(., alpha == 1))))))
qplot(x=x, y=y, color=fnct, data=df, geom="line",  log="x", ylim=c(.97, 1),
      main=bquote(pstable(x, alpha == 1.001, beta == 0.5)))+
  annotation_logticks(sides="b")
## but	alpha = 1   is only somewhat better as approximation:
showProc.time() #

c7 <- chk.pd.stable(1.2, -0.2,	 -40, 30)
c8 <- chk.pd.stable(1.5, -0.999, -40, 30)# two warnings

showProc.time() #

### Newly found -- Marius Hofert, 18.Sept. (via qstable):
stopifnot(all.equal(qstable(0.6, alpha = 0.5, beta = 1,
			    tol=1e-15),
			    2.63641788208582))    #This comes from qEstable
## was    2.636426573120147))
##--> which can be traced to the first of
stopifnot(pstable(q= -1.1, alpha=0.5, beta=1) == 0,
	  pstable(q= -2.1, alpha=0.6, beta=1) == 0)

## Found by Tobias Hudec, 2 May 2015:
stopifnot(
    all.equal(1.5, qstable(p=0.5, alpha=1.5, beta=0, gamma=2,   delta = 1.5)),
    all.equal(1.5, qstable(p=0.5, alpha=0.6, beta=0, gamma=0.2, delta = 1.5))
)

## Stable(alpha = 1/2, beta = 1, gamma, delta, pm = 1)	<===>  Levy(delta, gamma)
source(system.file("xtraR", "Levy.R", package = "stablecpp"), keep.source=interactive())
##-> dLevy(x, mu, c, log)                and
##-> pLevy(x, mu, c, log.p, lower.tail)

set.seed(101)
show.Acc <- (require("Rmpfr"))
if(show.Acc) { ## want to see accuracies, do not stop "quickly"
    format.relErr <- function(tt, cc)
        format(as.numeric(relErr(tt, cc)), digits = 4, width = 8)
}
##        pstable() is now as accurate as dstable()

pTOL <- 1e-15  # typically see relErr of  1 to 2 e-15
dTOL <- 1e-14 # typically see relErr of  1.3...3.9 e-15
showProc.time()

## Note that dstable() is more costly than pstable()
for(ii in 1:32) {
    Z <- rnorm(2)
    mu	<- sign(Z[1])*exp(Z[1])
    sc <- exp(Z[2])
    x <- seq(mu, mu+ sc* 100*rchisq(1, df=3),
	     length.out= 512)
    ## dLevy() and pLevy() using only pnorm() are "mpfr-aware":
    x. <- if(show.Acc) mpfr(x, 256) else x
    pL <- pLevy(x., mu, sc)
    pS <- pstable(x, alpha=1/2, beta=1, gamma=sc, delta=mu,
		  pm = 1)
    dL <- dLevy(x., mu, sc)
    dS <- dstable(x, alpha=1/2, beta=1, gamma=sc, delta=mu,
		  pm = 1)
    if(show.Acc) {
	cat("p: ", format.relErr(pL, pS), "\t")
	cat("d: ", format.relErr(dL, dS), "\n")
    } else {
	cat(ii %% 10)
    }
    stopifnot(all.equal(pL, pS, tol = pTOL),
	      all.equal(dL, dS, tol = dTOL))
}; cat("\n")

dev.off()

showProc.time()## ~ 75 sec  (doExtras=TRUE) on lynne (2012-09)
