require(stablecpp)
require(ggplot2)

##Functor which contains an approximation for stable p as a function of p for
#student t with df=alpha
pt_solve <- setClass("pt_solve",
                        slots = c(rplus = "numeric",
                                  rminus = "numeric",
                                  knot1 = "numeric",
                                  knot2 = "numeric",
                                  a0 = "numeric",
                                  a1 = "numeric",
                                  a2 = "numeric",
                                  a3 = "numeric",
                                  p = "numeric"
                                )
)

setMethod("initialize", "pt_solve",
          definition = function(.Object, pp, alpha, beta, lower_tail, log_p)
            {
              pi2<-pi/2
              c_stable_plus <- sin(pi2*alpha )*gamma(alpha)/pi*alpha*(1+beta)
              c_stable_minus <- sin(pi2*alpha )*gamma(alpha)/pi*alpha*(1-beta)
              c_t <- gamma((alpha+1)/2)/(sqrt(alpha*pi)*gamma(alpha/2))*alpha^((alpha+1)/2)
              .Object@rplus <- c_stable_plus/c_t;
              .Object@rminus <- c_stable_minus/c_t;
      # construct a cubic spline for the mapping of pt to cdf
             if (alpha<1 && beta==1){
               .Object@knot1 = pt(-tan(pi2*alpha), alpha)
             } else
               .Object@knot1 = .01
             if (alpha<1 && beta==-1) {
               .Object@knot2 = pt(tan(pi2*alpha), alpha)
             } else
               .Object@knot2 = .99
             dk=.Object@knot2-.Object@knot1
             .Object@a0 = .Object@rminus*.Object@knot1
             .Object@a1 = .Object@rminus
              b0 = 1-.Object@rplus*(1-.Object@knot2)-.Object@rminus*.Object@knot2
              b1 = .Object@rplus-.Object@rminus
             .Object@a2 = -(dk*b1-3*b0)/(dk*dk)
             .Object@a3 = (dk*b1-2*b0)/(dk*dk*dk)
             .Object@p = ifelse(log_p, exp(pp), pp)
             .Object@p = ifelse(lower_tail, .Object@p, 1-.Object@p)
             .Object
          }
    )

value <- function(pt_, ps) {pt_}

setMethod("value",
          signature = c(pt_ = "numeric",
                        ps = "pt_solve"),
          definition = function(pt_, ps) {
                         if (pt_<=ps@knot1)
                           ret <-pt_*ps@rminus-ps@p
                         else if (pt_>=ps@knot2)
                           ret <-1-ps@rplus*(1-pt_)-ps@p
                         else {
                           ptmk1 = pt_-ps@knot1
                           ret = (((ps@a3*ptmk1)+ps@a2)*ptmk1+ps@a1)*ptmk1+ps@a0 - ps@p
                         }
                         ret
                       }
         )

show_p_vs_pt <- function(alpha, beta) {
  pt_<-(1:499)/500
  q<-qt(pt_,alpha)
  p<-pstable(q, alpha, beta)
  df<-data.frame(pt_=pt_, q=q, p=p)
  df$type <- "exact"
  pt_s<-pt_solve(pp=0, alpha=alpha, beta=beta,
                 lower_tail=T, log_p=F)
  fit <- rep(NA, length(pt_))
  for (i in 1:length(pt_)) {
    fit[i] <- value(pt_[i], pt_s)
  }
  df<-rbind(df,data.frame(pt_=pt_, q=q, p=fit, type="approx"))
  gph<-qplot(x=pt_,y=p, data=df, geom="line",
             xlim=c(0,1), ylim=c(0,1), color=type)
  show(gph)
  df
}

q_guess <- function(p_vec, alpha, beta, lower_tail, log_p){
  ret <- rep(NA, length(p_vec))
  for (i in (1:length(p_vec))) {
    pt_s<-pt_solve(pp=p_vec[i], alpha=alpha, beta=beta,
                   lower_tail=lower_tail, log_p=log_p)
    r <- uniroot(value, interval=c(0,1), pt_s, tol = 1e-14, maxiter = 200);
    pt_ <- r$root
    ret[i] <- qt(max(0, min(pt_,1)), alpha)
  }
  ret
}

plot_guess<- function(alpha, beta, q_lim=c(-100,100)) {
  q <- seq(from=q_lim[1], to=q_lim[2], length.out=400)
  p0 <- pstable(q, alpha, beta)
  df<-data.frame(q=q, p0=p0, q_guess=q_guess(p0, alpha, beta, T, F))
  gph<-qplot(x=q, y=q_guess, data=df, color=I("black"), ylim=q_lim)+
    geom_line(aes(x=q, y=q), color="red")
  show(gph)
  df
}

q_guess_small_alpha <- function(p, alpha, beta, ylim=100) {
  pi2<-pi/2
  c_stable <- sin(pi2*alpha)*gamma(alpha)/pi
  p_low<-(1-beta)*.5*(1-exp(-1))*.75
  p_high<-1-(1+beta)*.5*(1-exp(-1))*.75
  p_mid <-.5 - .5*beta*(1-exp(-1))*.75
  q_low <- -(c_stable/(.5*(1-exp(-1))*.8))^(1/alpha)
  q_high <- -q_low
  q<-rep(NA,length(p))
  type<-rep(NA, length(p))
  sel <- p_low<p & p < p_high
  type[sel]<-"mid"
  adj<-3

  peak_stable <- gamma(1 + 1/alpha) /pi
  peak_student <- dt(0, alpha/adj)
  ratio <- peak_student/peak_stable
  q[sel]<-  qt(p[sel]+.5-p_mid, alpha/adj)*ratio
  q_low <- qt(p_low+.5-p_mid, alpha/adj)*ratio
  q_high <- qt(p_high+.5-p_mid, alpha/adj)*ratio
  sel<-p<p_low
  q[sel] <- -((c_stable*(1-beta))/p[sel])^(1/alpha)+q_low+((c_stable*(1-beta))/p_low)^(1/alpha)
  type[sel] <- "low"
  sel<-p>p_high
  q[sel] <- ((c_stable*(1+beta))/(1-p[sel]))^(1/alpha)+q_high-((c_stable*(1+beta))/(1-p_high))^(1/alpha)
  type[sel]<-"high"
  out <- data.frame(type=type, p=p, q=q)
  out$exact <- with(out,qstable(p, alpha, beta, pm=1))
  gph<-qplot(x=p,y=q, data=out, color=type, ylim=c(-ylim,ylim)) +
    geom_line(aes(x=p, y=exact), color = "black") +
    geom_vline(aes(xintercept=p_low), color="gray")+
    geom_vline(aes(xintercept=p_high),color="gray")
  show(gph)
  out
}

sigma_student<- function(alpha, method=1) {
  if (method==1){
    # rescale the student distribution so it's pdf at 0 matches the pdf of stable
    peak_stable <- gamma(1 + 1/alpha) /pi
    peak_student <- dt(0, alpha)
    ret <- peak_student/peak_stable
  } else {x
    # rescale the student distribution so it matches the tails of the cdf of
    # the stable distribution
    pi2 <- pi/2
    # the cdf tails of the standard distribution are asymptotic to c * x ^ (-alpha)
    c_student <<- gamma((alpha+1)/2)/(sqrt(alpha*pi)*gamma(alpha/2))*alpha^((alpha-1)/2)
    c_stable <<- sin(pi2*alpha )*gamma(alpha)/pi
    ret <- (c_student/c_stable)^(-1/alpha)
  }
  ret
  }

compare_pdf<-function(alpha, beta) {
  pi2 <- pi/2
  sigma_t <<- sigma_student(alpha)
  a <- -10*sigma_t
  b <- 10*sigma_t
  fact <<-(-a/.00001)^(1/100)
  x <<- (-100:100)/10*sigma_t
  n <- length(x)
  df <<- data.frame(type=rep("student",n),x = x, cdf=pt(x/sigma_t,alpha),
                    pdf=dt(x/sigma_t,alpha)/sigma_t )
  df <<- rbind(df,data.frame(type=rep("stable",n),x=x, cdf=pstable(x,alpha,beta),
                             pdf=dstable(x,alpha,beta)))
  gph_cdf <<-qplot(x=x, y=cdf, data=df, geom="line", color=type,
                   main=paste("Comparison of pstable(x,", format(alpha,digits=4),
                              ",", format(beta,digits=4),
                              ") to pt(x/sigma_t,", format(alpha,digits=4), ")"))
  gph_pdf <<-qplot(x=x, y=pdf, data=df, geom="line", color=type,
                   main=paste("Comparison of pstable(x,", format(alpha,digits=4),
                              ",", format(beta,digits=4),
                              ") to pt(x/sigma_t,", format(alpha,digits=4), ")"))
  df_p_vs_pt <-data.frame(pstable=df$cdf[df$type=="stable"], pt=df$cdf[df$type=="student"])
  gph_p_vs_pt <<- qplot(x=pt, y=pstable, data=df_p_vs_pt, geom="point",
                       main=paste("pstable(x,", format(alpha,digits=4),",",format(beta,digits=4),
                                  ") vs pstable(x/sigma_t,",format(alpha,digits=4),")"))
  show(gph_pdf)
  #show(gph_p_vs_pt)
}

student_asymptotic <- function(alpha, n, m, xlim=c(-10,10)) {
  # Display the distribution of large sums of variables w Student t distribution

  sigma_t <- sigma_student(alpha)

  # the scale of a sum is given by L(1/alpha) norm
  rtt<-sigma_t*rowSums(array(rt(n*m, alpha),dim=c(n,m)))/m^(1/alpha)
  rtt <- sort(rtt)
  r_ecdf <- ecdf(rtt)
  scaled_pt <- function(x, alpha) { pt(x/sigma_t, alpha)  }

  ks.stable <- ks.test(rtt, "pstable", alpha, 0)
  ks.student <- ks.test(rtt, "scaled_pt", alpha)

  df_note <-data.frame(x=rep(xlim[1]+.02*(xlim[2]-xlim[1]),2),y=c(1,.95, .85, .8),
                       label=c("Kolmogorov-Smirnov vs pstable",
                               paste("D =",format(ks.stable$statistic,digits=6),
                                     ", p-value =",format(ks.stable$p.value, digits=4)),
                                "Kolmogorov-Smirnov vs pt(x/sigma_t, alpha)",
                               paste("D =",format(ks.student$statistic,digits=6),
                                     ", p-value =",format(ks.student$p.value, digits=4))))

  df<-data.frame(x=rtt, empirical=r_ecdf(rtt))
  df_theory<-data.frame(type="stable",x=rtt, cdf=pstable(rtt, alpha, 0))
  df_theory<-rbind(df_theory, data.frame(type="student", x=rtt, cdf=pt(rtt/sigma_t,alpha)))
  gph <- qplot(x=x, y=empirical, data=df, color=I("blue"), xlim=xlim,
               ylab="cdf", main=paste("Distribution of",n, "Scaled Sums of rt(",m,", ", alpha,")"))+
          geom_line(aes(x=x, y=cdf, color=type),data=df_theory) +
          geom_rug(sides="b") +
          geom_label(aes(x=x,y=y,label=label),data=df_note, hjust="left")
  show(gph)
  stable_fit(rtt, type="mle")
}

show_tail_ratio <- function(alpha, beta, lower.tail=F) {
  c_stable <- sin(pi/2*alpha)*gamma(alpha)/pi
  q_low <- (c_stable/(.5*(1-exp(-1))))^(1/alpha)
  print(q_low)
  if (lower.tail) {
    c<-c_stable * (1 - beta)
  } else {
    c<-c_stable * (1+beta)
  }
  qmax <-(.Machine$double.xmax)^.5
  qmin <- q_low
  q_ratio <- qmin/qmax
  q <- qmax * q_ratio ^ ((0:200)/200)
  if (lower.tail)
    df<<-data.frame(x=q, stable=pstable(-q,alpha,beta,pm=1,lower.tail=T))
  else
    df <<-data.frame(x=q, stable=pstable(q,alpha,beta,pm=1,lower.tail=F))
  df$pareto<<- c * q^-alpha
  df$ratio<<-with(df, stable/pareto)
  gph<-qplot(x=pareto, y=ratio, data=df, geom="line")
  show(gph)
}
