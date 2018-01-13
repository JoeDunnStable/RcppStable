nolan_pm1_to_lindsey <- function(alpha, beta, gamma, delta) {
# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
tail <- alpha
loc <- delta
# the others require calcs:
del2 <- cos(pi/2 * alpha)^2 + (-beta)^2*sin(pi/2 * alpha)^2
del <- sqrt(del2) * sign(1-alpha)
eta_a <- min(alpha, 2-alpha)
# the lindsey-(3) beta:
skew <- 2/(pi*eta_a)*acos( cos(pi/2 * alpha) / del )
# the lindsey-(3) scale:
disp <- ( (del*gamma^alpha) / cos(pi/2 * alpha) )^(1/alpha)
data.frame(tail=tail, skew=skew, disp=disp, loc=loc)
}

pstable_lindsey<-function(q, alpha, beta, gamma, delta, eps=1e-6) {
  df<-nolan_pm1_to_lindsey(alpha, beta, gamma, delta)
  stable::pstable(q, tail = df$tail, skew=df$skew, disp =df$disp, loc  =df$loc, eps=eps)
}

dstable_lindsey<-function(q, alpha, beta, gamma, delta, eps=1e-6) {
  df<-nolan_pm1_to_lindsey(alpha, beta, gamma, delta)
  stable::dstable(q, tail = df$tail, skew=df$skew, disp =df$disp, loc  =df$loc, eps=eps)
}

qstable_lindsey<-function(p, alpha, beta, gamma, delta, eps=1e-6) {
  df<-nolan_pm1_to_lindsey(alpha, beta, gamma, delta)
  stable::qstable(p, tail = df$tail, skew=df$skew, disp =df$disp, loc  =df$loc, eps=eps)
}
