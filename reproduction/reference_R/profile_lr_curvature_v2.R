# Additive geometry and coupled bounds. Requires the unchanged v1 modules.
# Bounds below are justified in exact arithmetic. Base-R floating point is
# NOT outward-rounded interval arithmetic; never label it a machine certificate.

plr2_norm <- function(x) sqrt(sum(x*x))

plr2_three_fit <- function(x, counts, means) {
  stopifnot(length(x)==3L, length(counts)==3L, length(means)==3L,
            x[1]==0, all(diff(x)>0), all(counts>0), all(is.finite(c(x,counts,means))))
  b <- sqrt(counts/sum(counts))
  project <- function(v) v-b*sum(b*v)
  u <- x[2]/x[3]
  rays <- cbind(project(sqrt(counts)*c(0,u,1)),project(sqrt(counts)*c(0,1,1)))
  y <- project(sqrt(counts)*means)
  candidates <- rbind(c(0,0),c(max(0,sum(y*rays[,1])/sum(rays[,1]^2)),0),
                      c(0,max(0,sum(y*rays[,2])/sum(rays[,2]^2))))
  interior <- as.numeric(qr.solve(rays,y))
  if (all(interior>=0)) candidates <- rbind(candidates,interior)
  values <- apply(candidates,1,function(coef) sum((y-rays%*%coef)^2))
  i <- which.min(values)
  coef <- candidates[i,]
  delta <- sum(coef)
  ratio <- if (delta>0) (u*coef[1]+coef[2])/delta else 1
  t <- max(0,min(1,(1/ratio-1)/((1-u)/u)))
  a <- weighted.mean(means-delta*emax_shape(x,1-t,max(x)),counts)
  angle <- acos(max(-1,min(1,sum(rays[,1]*rays[,2])/
                               (plr2_norm(rays[,1])*plr2_norm(rays[,2])))))
  list(sse=values[i],a=a,delta=delta,t=t,angle=angle,
       coefficients=coef,analytic_active_set=TRUE,outward_rounding_certified=FALSE)
}

plr2_flat_tail <- function(cutoff,n,angle) {
  stopifnot(is.finite(cutoff),cutoff>=0,n>3,angle>0,angle<pi)
  k <- expm1(cutoff/n)
  nu <- n-3
  0.5*pf((nu+1)*k,1,nu+1,lower.tail=FALSE)+
    angle/(2*pi)*pf(nu*k/2,2,nu,lower.tail=FALSE)
}

plr2_three_full_error <- function(s,r) {
  stopifnot(s$m==3L,s$delta>0,r>0)
  Hsat <- plr2_vectors(s,0)$H
  Hlin <- plr2_vectors(s,1)$H
  orthogonal <- function(x,along) x-along*sum(x*along)/sum(along^2)
  nF <- orthogonal(Hlin,Hsat)
  nF <- nF/plr2_norm(nF)
  nOther <- orthogonal(Hsat,Hlin)
  nOther <- nOther/plr2_norm(nOther)
  nL <- orthogonal(s$w,s$v)
  zeta <- plr2_norm(nL)
  nL <- nL/zeta
  margins <- c(cone_linear=s$delta*sum(nOther*s$v),
               strip_amplitude=s$delta*plr2_norm(orthogonal(s$v,s$w)),
               strip_upper=(1-s$t)*zeta)
  eligible <- min(margins)>=r
  offsetF <- max(0,s$delta*sum(nF*s$v))
  offsetL <- s$t*zeta
  upper <- if (!eligible) r^2 else if (s$t==0) 0 else
    min(r^2,2*r*(r*plr2_norm(nF-nL)+abs(offsetF-offsetL)))
  list(eligible=eligible,upper=upper,margins=margins,
       same_halfspace=eligible && s$t==0,outward_rounding_certified=FALSE)
}

plr2_vectors <- function(s,t) {
  h <- emax_shape(s$x,1-t,max(s$x))
  hp <- -s$g/(1+s$g*t)^2
  list(H=s$project(sqrt(s$counts)*h/s$sigma),
       Hp=s$project(sqrt(s$counts)*hp/s$sigma))
}

plr2_shape_variation <- function(s,lo,hi) {
  mid <- (lo+hi)/2
  a <- (hi-lo)/2
  v <- plr2_vectors(s,mid)
  second <- plr2_norm(sqrt(s$counts)*2*s$g^2/(1+s$g*lo)^3)/s$sigma
  list(mid=mid,a=a,H=v$H,Hp=v$Hp,second=second,
       radius=a*plr2_norm(v$Hp)+a*a*second/2)
}

plr2_ball_interval <- function(direction,target,radius,limits) {
  # Coefficients d for ||d*direction-target|| <= radius, intersected with limits.
  den <- sum(direction^2)
  center <- sum(direction*target)/den
  residual <- target-center*direction
  slack <- radius^2-sum(residual^2)
  cushion <- 64*.Machine$double.eps*max(1,radius^2,sum(target^2))
  if (slack < -cushion) return(NULL)
  half <- sqrt(max(0,slack+cushion)/den)
  ab <- c(max(limits[1],center-half),min(limits[2],center+half))
  if (ab[1]>ab[2]) NULL else ab
}

plr2_delta_limits <- function(s,r) {
  d <- 2*r*s$sigma/sqrt(s$n)*sqrt(1/s$allocation[1]+1/tail(s$allocation,1))
  c(max(0,s$delta-d),s$delta+d)
}

plr2_search <- function(enclose,lo=0,hi=1,tol=1e-3,max_intervals=256L,
                        maximize=FALSE) {
  stopifnot(tol>0,max_intervals>=1L)
  first <- enclose(lo,hi)
  boxes <- matrix(c(lo,hi,first$bound),nrow=1)
  best <- first$value
  arg <- first$arg
  evaluations <- 1L
  repeat {
    j <- if (maximize) which.max(boxes[,3]) else which.min(boxes[,3])
    outer <- boxes[j,3]
    gap <- if (maximize) outer-best else best-outer
    if (gap<=tol || nrow(boxes)>=max_intervals) break
    ab <- boxes[j,1:2]
    mid <- mean(ab)
    if (mid==ab[1] || mid==ab[2]) break
    left <- enclose(ab[1],mid)
    right <- enclose(mid,ab[2])
    for (candidate in list(left,right)) {
      better <- if (maximize) candidate$value>best else candidate$value<best
      if (better) { best <- candidate$value; arg <- candidate$arg }
    }
    boxes <- rbind(boxes[-j,,drop=FALSE],c(ab[1],mid,left$bound),
                   c(mid,ab[2],right$bound))
    evaluations <- evaluations+2L
  }
  outer <- if (maximize) max(boxes[,3]) else min(boxes[,3])
  lower <- if (maximize) best else min(outer,best)
  upper <- if (maximize) max(outer,best) else best
  list(lower=lower,upper=upper,arg=arg,gap=upper-lower,
       converged=upper-lower<=tol,intervals=nrow(boxes),evaluations=evaluations,
       bound_in_exact_arithmetic=TRUE,outward_rounding_certified=FALSE)
}

plr2_coupled_remainder <- function(s,r,null=FALSE,tol=1e-3,max_intervals=256L) {
  limits <- plr2_delta_limits(s,r)
  tau <- s$delta*(1-plr_deficit_v1(s$t))
  delta_at <- function(t) tau/(1-plr_deficit_v1(t))
  remainder <- function(delta,t,H) delta*(H-s$v)-s$delta*(t-s$t)*
    plr2_vectors(s,s$t)$Hp
  enclose <- function(lo,hi) {
    b <- plr2_shape_variation(s,lo,hi)
    local_limits <- limits
    if (null) local_limits <- c(max(limits[1],delta_at(lo)),min(limits[2],delta_at(hi)))
    if (local_limits[1]>local_limits[2])
      return(list(bound=0,value=0,arg=c(delta=s$delta,t=s$t)))
    full <- plr2_ball_interval(b$H,s$delta*s$v,2*r+limits[2]*b$radius,local_limits)
    linear <- plr2_ball_interval(s$v,s$delta*s$v-(b$mid-s$t)*s$w,
                                2*r+b$a*plr2_norm(s$w),local_limits)
    if (is.null(full) && is.null(linear))
      return(list(bound=0,value=0,arg=c(delta=s$delta,t=s$t)))
    ends <- unlist(list(full,linear),use.names=FALSE)
    ends <- range(ends)
    upper <- max(vapply(ends,function(d) {
      plr2_norm(remainder(d,b$mid,b$H))+b$a*plr2_norm(d*b$Hp-s$w)
    },numeric(1)))+ends[2]*b$second*b$a^2/2
    # Feasible center/endpoints give a lower bound on the supremum.
    best <- 0
    arg <- c(delta=s$delta,t=s$t)
    for (t in unique(c(lo,b$mid,hi,s$t[s$t>=lo & s$t<=hi]))) {
      H <- plr2_vectors(s,t)$H
      ab1 <- plr2_ball_interval(H,s$delta*s$v,2*r,limits)
      ab2 <- plr2_ball_interval(s$v,s$delta*s$v-(t-s$t)*s$w,2*r,limits)
      ds <- if (null) delta_at(t) else unique(unlist(lapply(list(ab1,ab2),function(ab) {
        if (is.null(ab)) return(NULL)
        c(mean(ab),ab+(mean(ab)-ab)*1e-10)
      })))
      for (d in ds) {
        f <- d*H-s$delta*s$v
        ell <- (d-s$delta)*s$v+(t-s$t)*s$w
        if (min(plr2_norm(f),plr2_norm(ell))<=2*r) {
          val <- plr2_norm(f-ell)
          if (val>best) { best <- val; arg <- c(delta=d,t=t) }
        }
      }
    }
    list(bound=max(upper,best),value=best,arg=arg)
  }
  plr2_search(enclose,tol=tol,max_intervals=max_intervals,maximize=TRUE)
}

plr2_profile_bounds <- function(s,z,null=TRUE,linear=FALSE,tol=1e-3,
                               max_intervals=256L,extra_t=numeric()) {
  z <- s$project(z)
  r <- plr2_norm(z)
  stopifnot(r>0)
  limits <- plr2_delta_limits(s,r)
  tau <- s$delta*(1-plr_deficit_v1(s$t))
  delta_at <- function(t) tau/(1-plr_deficit_v1(t))
  at <- function(t) {
    H <- plr2_vectors(s,t)$H
    d <- if (null) delta_at(t) else max(0,sum(H*(z+s$delta*s$v))/sum(H^2))
    f <- if (linear) (d-s$delta)*s$v+(t-s$t)*s$w else d*H-s$delta*s$v
    list(value=sum((z-f)^2),arg=c(delta=d,t=t))
  }
  enclose <- function(lo,hi) {
    b <- plr2_shape_variation(s,lo,hi)
    if (null) {
      ds <- c(delta_at(lo),delta_at(hi))
      dc <- mean(ds)
      if (linear) {
        center <- (dc-s$delta)*s$v+(b$mid-s$t)*s$w
        radius <- diff(ds)*plr2_norm(s$v)/2+b$a*plr2_norm(s$w)
      } else {
        center <- dc*b$H-s$delta*s$v
        radius <- diff(ds)*plr2_norm(b$H)/2+ds[2]*b$radius
      }
      lower <- max(0,plr2_norm(z-center)-radius)^2
    } else {
      stopifnot(!linear)
      d <- max(limits[1],min(limits[2],sum(b$H*(z+s$delta*s$v))/sum(b$H^2)))
      distance <- plr2_norm(z+s$delta*s$v-d*b$H)
      lower <- max(0,distance-limits[2]*b$radius)^2
    }
    trials <- lapply(unique(c(lo,b$mid,hi,extra_t[extra_t>=lo & extra_t<=hi],
                               s$t[s$t>=lo & s$t<=hi])),at)
    best <- trials[[which.min(vapply(trials,`[[`,numeric(1),"value"))]]
    list(bound=min(lower,best$value),value=best$value,arg=best$arg)
  }
  plr2_search(enclose,tol=tol,max_intervals=max_intervals)
}

plr2_loss_increment <- function(s,z,theta) {
  H <- plr2_vectors(s,theta["t"])$H
  f <- theta["delta"]*H-s$delta*s$v
  ell <- (theta["delta"]-s$delta)*s$v+(theta["t"]-s$t)*s$w
  R <- f-ell
  -2*sum((z-ell)*R)+sum(R^2)
}

plr2_signed_sandwich <- function(s,z,A0,A1,H0,H1) {
  D <- function(fit) plr2_loss_increment(s,z,fit$arg)
  lower <- D(A0)-D(H1)-A0$gap-H1$gap
  upper <- D(H0)-D(A1)+H0$gap+A1$gap
  list(lower=lower,upper=upper,width=upper-lower,
       bound_in_exact_arithmetic=TRUE,outward_rounding_certified=FALSE)
}
