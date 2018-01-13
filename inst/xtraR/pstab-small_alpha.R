require("stablecpp")
require("ggplot2")
require("plyr")

pPareto <- stablecpp:::pPareto

pdf("pstab-small_alpha.pdf")

## Now *really* small alpha --> 0:
##     --------------------------

pstable_alpha0 <-function(x, beta, lower.tail=T, log_p=F) {
  ret<-rep(NA, length(x))
  ret <- ifelse(x<0,.5*(1-beta)*(1-exp(-1)),
                    1-.5*(1+beta)*(1-exp(-1)))
  ret <- ifelse(rep(lower.tail,length(x)), ret, 1-ret)
  ret <- ifelse(rep(log_p, length(x)), log(ret), ret)
  ret
}

pdstab.alpha <- function(x, betas, alphas = 2^-(40:10),title, ...)
{
  stopifnot(is.numeric(x),
            is.numeric(betas), -1 <= betas, betas <= 1,
            is.numeric(alphas), 0 <= alphas, alphas <= 2)
  if(is.unsorted(alphas)) alphas <- sort(alphas)
  if(is.unsorted(beta)) betas<-sort(betas)
  one_ab<-function(df){
    data.frame(x=x,cdf=pstable(x, alpha=df$alpha, beta = df$beta, pm = 1),
               cdf0=pstable_alpha0(x, df$beta))
  }
  df<-expand.grid(alpha=alphas,beta=betas)
  df_out<<-ddply(df,.(alpha,beta),one_ab)
  qplot(x=alpha,y=abs((cdf)-(cdf0)),data=df_out,color=as.factor(beta),
        linetype=as.factor(x), geom="line",log="xy",
        xlab = quote(alpha),
        ylab = expression(abs(F(x, alpha)-F(x, alpha=0))),
        main = title)+
    scale_color_discrete(guide_legend(title="beta"))+
    scale_linetype_discrete(guide_legend(title="x"))
}

## vary beta,  keep x :
pdstab.alpha(x = 10, betas = c(0.1,.3,.5,.7,.9),
             title=expression(pstable(x == -1, alpha, beta, pm == 0)-pstable(x == -1, alpha == 0, beta)))
## vary x,  keep beta
pdstab.alpha(x =  1:10, beta = 0.3,
             title=expression(pstable(x, beta == .3, alpha, pm == 0)-pstable(x, beta == .3, alpha == 0, pm == 0)))

## Plots suggest a simple formula
##  log(f(x, alpha, beta))=  c(x,beta) + log(alpha)
##      f(x, alpha, beta) =  C(x,beta) * (alpha)   -- ???

## for normal alpha, it looks different {which is reassuring!}
pdstab.alpha(x = -1, betas = 0.3, alphas = seq(1/128, 2, length=100),
             title=expression(pstable(x == -1, alpha, beta, pm == 0)-pstable(x == -1, 0, beta, pm == 0)))

dev.off()
stop()

## Show the formula, we derived:
## f(x, \alpha, \beta) ~  \alpha * \frac{1 + \beta}{2e \abs{x + \pi/2 \alpha\beta}}
## ONLY ok, when  x > zeta := - beta * tan(pi/2 *alpha)
## otherwise, "swap sign" of (x, beta)
cdf.smlA <- function(x, alpha, beta, log = FALSE) {
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

alpha <- 1e-15 ## must be larger than cutoff in pstable
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

dev.off()
