library(RColorBrewer)

#############################################################
# Functions used to calculate the theoretical reliability diagrams (Figure 3)
library(GoFKernel) # inverse()

# Forecast parameters
# unfocused
c_unf = 1.5
# lopsided
c_lop = 0.7
# piecewise uniform
c_pwUnif = 0.5

# Conditional non exceedance probabilities
# unfocused
cond_nep_unf = function(alpha, t = 0, c = c_unf){
  Phi_plus = function(x) 0.5*(pnorm(x) + pnorm(x-c))
  Phi_minus = function(x) 0.5*(pnorm(x) + pnorm(x+c))
  inverse_Phi_plus = inverse(Phi_plus)
  inverse_Phi_minus = function(u) inverse_Phi_plus(u) - c
  
  weight_plus = dnorm(t - inverse_Phi_plus(alpha))
  weight_minus = dnorm(t - inverse_Phi_minus(alpha))
  
  return((weight_plus*pnorm(inverse_Phi_plus(alpha)) + weight_minus*pnorm(inverse_Phi_minus(alpha)))/(weight_plus + weight_minus))
}

# lopsided
cond_nep_lop = function(alpha, t = 0, s = c_lop){
  mu = function(alpha,delta){
    ifelse(delta == 1, t - qnorm((alpha+delta)/(1+delta)),
           ifelse(delta == -1, t - qnorm(alpha/(1-delta)),
                  ifelse(alpha <= (1-delta)/2, t - qnorm(alpha/(1-delta)), t - qnorm((alpha+delta)/(1+delta)))))
  }
  weight_plus = dnorm(mu(alpha,s))
  weight_minus = dnorm(mu(alpha,-s))
  
  return((weight_plus*pnorm(t-mu(alpha,s)) + weight_minus*pnorm(t-mu(alpha,-s)))/(weight_plus + weight_minus))
}

# piecewise uniform
cond_nep_pwUnif = function(alpha, t = 0, sigma = c_pwUnif){
  inverse_G1 = function(u) if(u <= 1/2) 2*u else 4*(u-1/2) + 1
  inverse_G2 = function(u) if(u <= 1/4) 4*u else if(u <= 3/4) 2*(u - 1/4)+1 else 4*(u-3/4) + 2
  inverse_G3 = function(u) if(u <= 1/2) 4*u else 2*(u - 1/2)+2
  
  weight1 = dnorm((t - inverse_G1(alpha))/sigma)
  weight2 = dnorm((t - inverse_G2(alpha))/sigma)
  weight3 = dnorm((t - inverse_G3(alpha))/sigma)
  
  ptrue1 = function(x) if(x <= 0) 0 else if(x <= 1) x/2 else if(x <= 2) 1/2 + 0.1*(x-1) else if(x <= 3) 3/5 + 2/5*(x-2) else 1
  ptrue2 = function(x) if(x <= 0) 0 else if(x <= 1) x/10 else if(x <= 2) 1/10 + 0.8*(x-1) else if(x <= 3) 9/10 + 1/10*(x-2) else 1
  ptrue3 = function(x) if(x <= 0) 0 else if(x <= 1) 2*x/5 else if(x <= 2) 2/5 + 0.1*(x-1) else if(x <= 3) 1/2 + 1/2*(x-2) else 1
  
  return((weight1*ptrue1(inverse_G1(alpha)) + weight2*ptrue2(inverse_G2(alpha)) + weight3*ptrue3(inverse_G3(alpha)))/(weight1 + weight2 + weight3))
}

# Quantiles
# unfocused
cond_quant_unf = function(q,alpha = 0.5,c = c_unf){
  Phi_plus = function(x) 0.5*(pnorm(x) + pnorm(x-c))
  Phi_minus = function(x) 0.5*(pnorm(x) + pnorm(x+c))
  inverse_Phi_plus = inverse(Phi_plus)
  inverse_Phi_minus = function(u) inverse_Phi_plus(u) - c
  
  weight_plus = dnorm(q - inverse_Phi_plus(alpha))
  weight_minus = dnorm(q - inverse_Phi_minus(alpha))
  
  cond_cdf = function(t) (weight_plus*pnorm(t - q + inverse_Phi_plus(alpha)) + weight_minus*pnorm(t - q + inverse_Phi_minus(alpha)))/(weight_plus + weight_minus)
  
  # Estimate the quantile (inverse() to slow here...)
  upp = ceiling(qnorm(alpha,mean = q - inverse_Phi_minus(alpha)))
  low = floor(qnorm(alpha,mean = q - inverse_Phi_plus(alpha)))
  for(i in 1:5){
    ind = which(cond_cdf(seq(low,upp,10^-i))> alpha)[1]
    if(is.na(ind)) return(upp)
    upp = seq(low,upp,10^-i)[ind]
    low = seq(low,upp,10^-i)[ind-1]
  }
  return(low)
}

# lopsided
cond_quant_lop = function(q,alpha = 0.5,s = c_lop){
  mu = function(alpha,delta){
    ifelse(delta == 1, q - qnorm((alpha+delta)/(1+delta)),
           ifelse(delta == -1, q - qnorm(alpha/(1-delta)),
                  ifelse(alpha <= (1-delta)/2, q - qnorm(alpha/(1-delta)), q - qnorm((alpha+delta)/(1+delta)))))
  }
  weight_plus = dnorm(mu(alpha,s))
  weight_minus = dnorm(mu(alpha,-s))
  
  cond_cdf = function(t) (weight_plus*pnorm(t-mu(alpha,s)) + weight_minus*pnorm(t-mu(alpha,-s)))/(weight_plus + weight_minus)
  
  # Estimate the quantile (inverse() to slow here...)
  low = floor(min(qnorm(alpha,mean = mu(alpha,s)),qnorm(alpha,mean = mu(alpha,-s))))
  upp = ceiling(max(qnorm(alpha,mean = mu(alpha,s)),qnorm(alpha,mean = mu(alpha,-s))))
  for(i in 1:5){
    ind = which(cond_cdf(seq(low,upp,10^-i))> alpha)[1]
    if(is.na(ind) || low == upp) return(upp)
    upp = seq(low,upp,10^-i)[ind]
    low = seq(low,upp,10^-i)[ind-1]
  }
  return(low)
}

# piecewise uniform
cond_quant_pwUnif = function(q,alpha = 0.5,sigma = c_pwUnif){
  inverse_G1 = function(u) if(u <= 1/2) 2*u else 4*(u-1/2) + 1
  inverse_G2 = function(u) if(u <= 1/4) 4*u else if(u <= 3/4) 2*(u - 1/4)+1 else 4*(u-3/4) + 2
  inverse_G3 = function(u) if(u <= 1/2) 4*u else 2*(u - 1/2)+2
  
  weight1 = dnorm((q - inverse_G1(alpha))/sigma)
  weight2 = dnorm((q - inverse_G2(alpha))/sigma)
  weight3 = dnorm((q - inverse_G3(alpha))/sigma)
  
  ptrue1 = function(x) ifelse(x <= 0, 0,
                              ifelse(x <= 1, x/2,
                                     ifelse(x <= 2, 1/2 + 0.1*(x-1),
                                            ifelse(x <= 3, 3/5 + 2/5*(x-2), 1))))
  ptrue2 = function(x) ifelse(x <= 0, 0,
                              ifelse(x <= 1, x/10,
                                     ifelse(x <= 2, 1/10 + 0.8*(x-1),
                                            ifelse(x <= 3, 9/10 + 1/10*(x-2), 1))))
  ptrue3 = function(x) ifelse(x <= 0, 0, 
                              ifelse(x <= 1, 2*x/5,
                                     ifelse(x <= 2, 2/5 + 0.1*(x-1),
                                            ifelse(x <= 3, 1/2 + 1/2*(x-2), 1))))
  
  cond_cdf = function(t) (weight1*ptrue1(t-q+inverse_G1(alpha)) + weight2*ptrue2(t-q+inverse_G2(alpha)) + weight3*ptrue3(t-q+inverse_G3(alpha)))/(weight1 + weight2 + weight3)
  
  # Estimate the quantile (inverse() to slow here...)
  low = q-3
  upp = q+3
  for(i in 1:5){
    ind = which(cond_cdf(seq(low,upp,10^-i))> alpha)[1]
    if(is.na(ind)) return(upp)
    upp = seq(low,upp,10^-i)[ind]
    low = seq(low,upp,10^-i)[ind-1]
  }
  return(low)
}

# Moments
# unfocused
cond_moms_unf = list(
  cond_mom1 = function(m, c = c_unf){
    weight_plus = dnorm(m + c/2)
    weight_minus = dnorm(m - c/2)
    return((weight_plus*(m + c/2) + weight_minus*(m - c/2))/(weight_plus + weight_minus))
  },
  cond_mom2 = function(m, c = c_unf){
    eta = c(1,-1,1,-1)*c
    pm = c(1,1,-1,-1)
    mu = -eta/2+pm*sqrt(m - 1 - eta^2/4)
    weights = dnorm(mu)
    moments = mu^2 + 1
    return(sum(weights*moments)/sum(weights))
  },
  cond_mom3 = function(m, c = c_unf, eps = 10^-10){
    roots = matrix(NA,nrow = 2,ncol = 3)
    eta = c(c,-c)
    for(i in 1:2) roots[i,] = polyroot(c(0.5*(eta[i]^3 + 3*eta[i]) - m, 3*(eta[i]^2/2 + 1), 1.5*eta[i], 1))
    mu = as.vector(Re(roots[abs(Im(roots)) < eps])) # the real-valued roots
    weights = dnorm(mu)
    moments = mu^3 + 3*mu
    return(sum(weights*moments)/sum(weights))
  }
)

# lopsided
cond_moms_lop = list(
  cond_mom1 = function(m, s = c_lop){
    weight_plus = dnorm(m + 2*s*dnorm(0))
    weight_minus = dnorm(m - 2*s*dnorm(0))
    return((weight_plus*(m + 2*s*dnorm(0)) + weight_minus*(m - 2*s*dnorm(0)))/(weight_plus + weight_minus))
  },
  cond_mom2 = function(m, s = c_lop){
    delta = c(1,-1,1,-1)*s
    pm = c(1,1,-1,-1)
    mu = -2*delta*dnorm(0) + pm * sqrt(4*delta^2 *dnorm(0)^2 + m - 1)
    weights = dnorm(mu)
    moments = mu^2 + 1
    return(sum(weights*moments)/sum(weights))
  },
  cond_mom3 = function(m, s = c_lop, eps = 10^-10){
    roots = matrix(NA,nrow = 2,ncol = 3)
    delta = c(s,-s)
    for(i in 1:2) roots[i,] = polyroot(c(4*delta[i]*dnorm(0) - m, 3, 6*delta[i]*dnorm(0), 1))
    mu = as.vector(Re(roots[abs(Im(roots)) < eps])) # the real-valued roots
    weights = dnorm(mu)
    moments = mu^3 + 3*mu
    return(sum(weights*moments)/sum(weights))
  }
)

# piecewise uniform
cond_moms_pwUnif = list(
  cond_mom1 = function(m, sigma = c_pwUnif){
    mean_G = function(i) 1 + i/4
    mean_tildeG = function(i) ifelse(i == 1, 14,
                                     ifelse(i == 2, 15,
                                            ifelse(i == 3, 16,NA)))/10
    
    weight = function(i) dnorm((m - mean_G(i))/sigma)
    
    return(sum(weight(1:3)*(mean_tildeG(1:3) + (m-mean_G(1:3))))/sum(weight(1:3)))
  },
  cond_mom2 = function(m, sigma = c_pwUnif){
    p = (matrix(1,nrow = 3,ncol = 3) + diag(3))/4
    q = matrix(c(5,1,4,1,8,1,4,1,5),nrow = 3)/10
    
    a = function(p) colSums(p*(2*1:3 - 1))
    b = function(p) colSums(p*(3*(1:3)^2 - 3*1:3 + 1))/3
    
    pm = c(1,1,1,-1,-1,-1)
    
    mu = - rep(a(p),2)/2 + pm * sqrt((rep(a(p),2)/2)^2 - rep(b(p),2) + m)
    
    weights = dnorm(mu/sigma)
    
    moments = mu^2 + mu * a(q) + b(q)
    return(sum(weights*moments)/sum(weights))
  },
  cond_mom3 = function(m, sigma = c_pwUnif, eps = 10^-10){
    p = (matrix(1,nrow = 3,ncol = 3) + diag(3))/4
    q = matrix(c(5,1,4,1,8,1,4,1,5),nrow = 3)/10
    
    a = function(p) 3/2*colSums(as.matrix(p*(2*1:3 - 1)))
    b = function(p) colSums(as.matrix(p*(3*(1:3)^2 - 3*1:3 + 1)))
    d = function(p) colSums(as.matrix(p*((1:3)^4 - (1:3 - 1)^4)))/4
    
    mu = rep(NA,3)
    for(i in 1:3){
      roots = polyroot(c(d(p[,i]) - m, b(p[,i]), a(p[,i]), 1))
      mu[i] = Re(roots[abs(Im(roots)) < eps]) # the real-valued root
    }
    weights = dnorm(mu/sigma)
    moments = mu^3 + a(q)*mu^2 + b(q)*mu + d(q)
    return(sum(weights*moments)/sum(weights))
  }
)

# Plotting routines
plot_thresh = function(cond_nep,ts,main = "",xlab = expression(x == F(t)),ylab = expression(x[rc]), colors){
  cols = colors[c(2,4,5,6,8)]
  ltys = c(2,2,1,2,2)
  leg1 = 5:4
  leg2 = 3:1
  u = seq(0,1,length.out = 200)
  
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main = main,xlab = xlab,ylab = ylab)
  abline(a = 0,b = 1,lty = 1,col = "grey")
  
  for(i in 1:length(ts)){
    nep_values = sapply(u,cond_nep,t = ts[i])
    points(u,nep_values,type = 'l',lty = ltys[i],col = cols[i])
  }
  legend("topleft", legend = ts[leg1],lty = ltys[leg1],col = cols[leg1],title = expression(t), bg = "white")
  legend("bottomright", legend = ts[leg2],lty = ltys[leg2],col = cols[leg2], bg = "white")
}

plot_quant = function(cond_quant,alphas,main = NULL,xlab = expression(x = q[alpha](F)),ylab = expression(x[rc]),lim = c(-3,3), colors){
  cols = colors
  ltys = c(2,2,2,2,1,2,2,2,2)
  leg1 = 9:6
  leg2 = 5:1
  q = seq(lim[1],lim[2],length.out = 200)
  
  plot(NULL,xlim = lim,ylim = lim,main = main,xlab = xlab,ylab = ylab)
  abline(a = 0,b = 1,lty = 1,col = "grey")
  
  for(i in 1:length(alphas)){
    quant_values = sapply(q,cond_quant,alpha = alphas[i])
    points(q,quant_values,type = 'l',lty = ltys[i],col = cols[i])
  }
  legend("topleft", legend = alphas[leg1],lty = ltys[leg1],col = cols[leg1],title = expression(alpha))
  legend("bottomright", legend = alphas[leg2],lty = ltys[leg2],col = cols[leg2])
}

plot_moms = function(cond_moms,param = NULL,main = NULL,xlab = NULL,ylab = NULL,
                     lim = c(-5,5), colors, root_transform = TRUE){
  cols = colors[c(5,2,8)]
  ltys = c(1,2,2)
  leg = 3:1
  ns = 1:3
  
  if(is.null(xlab)){
    if(root_transform) xlab = expression(x == sqrt(m[n](F),n))
    else xlab = expression(x == m[n](F))
  }
  if(is.null(ylab)) ylab = expression(x[rc])
  
  cbrt = function(x) ifelse(x < 0, -(-x)^(1/3), x^(1/3))
  root = function(x,n) if(n == 1 || !root_transform) x else if(n == 2) sqrt(x) else cbrt(x)
  pow = function(x,n) if(n == 1 || !root_transform) x else x^n
  m = seq(lim[1],lim[2],length.out = 200)
  
  plot(NULL,xlim = lim,ylim = lim,main = main,xlab = xlab,ylab = ylab)
  abline(a = 0,b = 1,lty = 1,col = "grey")
  
  mom1_values = sapply(m,cond_moms[[1]])
  points(m,mom1_values,type = 'l',col = cols[1],lty = ltys[1])
  
  mom2_values = sapply(pow(m[m>0],2),cond_moms[[2]])
  points(m[m>0],root(mom2_values,2),type = 'l',col = cols[2],lty = ltys[2])
  
  mom3_values = sapply(pow(m,3),cond_moms[[3]])
  points(m,root(mom3_values,3),type = 'l',col = cols[3],lty = ltys[3])
  
  legend("bottomright", legend = ns[leg],lty = ltys[leg],col = cols[leg],title = expression(n))
}

#############################################################
# Empirical Reliability Diagrams
# Functions to plot various reliability diagrams
library(isotone)

# Marginal Reliability Diagram
# fcast needs to be a list with the following elements:
# fcast$F: a function, which maps a vector x of the same length as y to the forecast CDF values at the values x,
#   i.e., fcast$F(x)[i] = F_i(x_i)
# fcast$sample: a function with an argument cases, which samples from the forecast CDFs specified in cases, 
#   e.g., for cases = 1:n the function fcast$sample should sample once from each forecast
margreldiag = function(fcast,y,resampling = TRUE,n_resamples = 1000,region_level = 0.9){
  n_pts = 1000 # the max number of points used to determine the curve
  y = sort(y)
  if(length(y) > n_pts) t = sort(c(head(y,1),sample(y,n_pts),tail(y,1)))
  else t = y
  
  avg_fcast_y = sapply(t,function(x) mean(fcast$F(x)))
  marg_dist_y = ecdf(y)(t)
  
  plot(NULL,xlim = 0:1,ylim = 0:1,main = "",ylab = "sample NEP",xlab = "average forecast NEP")

  if(resampling){
    low = floor(n_resamples * (1-region_level)/2)
    up = n_resamples - low
    
    # Resampling from the average forecast distribution is the same as sampling from randomly drawn forecasts (cases)
    resampled_cases = matrix(sample(1:length(y),size = n_resamples*length(y),replace = TRUE),ncol = n_resamples)
    resamples = apply(resampled_cases,2,function(cases) fcast$sample(cases))
    
    marg_dist_y_resamples = apply(resamples, 2, function(s) ecdf(s)(t))
    marg_dist_y_resamples_sorted = apply(marg_dist_y_resamples,1,sort) # sorted pointwise
    
    polygon(c(avg_fcast_y,rev(avg_fcast_y)),
            c(marg_dist_y_resamples_sorted[low,],rev(marg_dist_y_resamples_sorted[up,])),
            border = NA,col = "lightblue1")
    points(avg_fcast_y,marg_dist_y_resamples_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(avg_fcast_y,marg_dist_y_resamples_sorted[up,],type = "l",lty = 1,col = "lightblue2")
  }
  
    abline(a=0,b=1,lty=2,col="grey")
    points(avg_fcast_y,marg_dist_y,type = "l")
}

# PIT-Reliability Diagram
# fcast needs to be a list with the following element:
# fcast$F: a function, which maps a vector x of the same length as y to the forecast CDF values at the values x,
#   i.e., fcast$F(x)[i] = F_i(x_i)
PITreldiag = function(fcast,y,resampling = TRUE,n_resamples = 1000,region_level = 0.9){
  # auxiliary function
  color_step_poly = function(lim,x_up,x_low,y_up,y_low,col = "lightblue1"){
    polygon(c(lim[1],rbind(x_up,c(x_up[-1],lim[2])),lim[2],rbind(rev(x_low),rev(x_low)),lim[1]),
            c(lim[1],c(rbind(y_up,y_up)),lim[2],rbind(rev(y_low),c(rev(y_low)[-1],lim[1])),lim[1]),
            border = NA,col = "lightblue1")
  }
  
  dist_PIT = ecdf(fcast$F(y))
  
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main = "", 
       xlab = expression(z),ylab = expression("fraction of PIT-values" <= z))
  
  if(resampling){
    low = floor(n_resamples * (1-region_level)/2)
    up = n_resamples - low

    resamples = sapply(1:n_resamples,function(i) runif(length(y)))
    
    t = seq(0,1,0.001)
    dist_resamples_t = apply(resamples, 2, function(s) ecdf(s)(t))
    dist_resamples_t_sorted = apply(dist_resamples_t,1,sort)
    
    color_step_poly(c(0,1),t,t,dist_resamples_t_sorted[up,],dist_resamples_t_sorted[low,])
    # polygon(c(0,t,1,rev(t),0),
    #         c(0,dist_resamples_t_sorted[low,],1,rev(dist_resamples_t_sorted[up,]),0),border = NA,col = "lightblue1")
    points(t, dist_resamples_t_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(t, dist_resamples_t_sorted[up,],type = "l",lty = 1,col = "lightblue2")
  }
  
  points(c(0,1),c(0,1),type = "l",lty = 2,col = "grey")
  plot(dist_PIT,do.points = FALSE,xlim = c(0,1),col.01line = NULL,verticals = TRUE,add=TRUE)
}

# Threshold Reliability Diagram
# fcast needs to be a list with the following element:
# fcast$F: a function, which maps a vector x of the same length as y to the forecast CDF values at the values x,
#   i.e., fcast$F(x)[i] = F_i(x_i)
threshreldiag = function(fcast,y,t,resampling = TRUE,n_resamples = 1000,region_level = 0.9,
                         digits = 3,inset_hist = TRUE,main = "",xlab = NULL,ylab = NULL, ...){
  if(is.null(xlab)) xlab = bquote(x == F(.(t)))
  if(is.null(ylab)) ylab = expression(x[rc])
  
  # binarize problem
  x = fcast$F(t)
  y = as.numeric(y <= t)
  sample = function() as.numeric(runif(y) <= x)
  score = function(x,y) mean((x-y)^2)
    
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main = main,xlab = xlab,ylab = ylab)
  
  x_rc = isoreg(x,y)$yf # values correspond to ORDERED forecast values!
  
  # compute score decomposition
  s = score(x,y)
  s_rc = score(x_rc,y[order(x)])
  s_mg = score(mean(y),y)
  mcb = s - s_rc
  dsc = s_mg - s_rc
  unc = s_mg
  
  if(resampling){
    low = floor(n_resamples * (1-region_level)/2)
    up = n_resamples - low
    pval_digits = ceiling(log(n_resamples,10))
    
    resamples = sapply(1:n_resamples,function(i) sample())
    x_rc_resamples = apply(resamples, 2, function(y) isoreg(x,y)$yf)
    x_rc_resamples_sorted = apply(x_rc_resamples,1,sort)
    
    ran_x = range(x)
    polygon(c(ran_x[1],sort(x),ran_x[2],rev(sort(x)),ran_x[1]),
            c(ran_x[1],x_rc_resamples_sorted[up,],ran_x[2],rev(x_rc_resamples_sorted[low,]),ran_x[1]),
            border = NA,col = "lightblue1")
    points(sort(x),x_rc_resamples_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(sort(x),x_rc_resamples_sorted[up,],type = "l",lty = 1,col = "lightblue2")
    box()
    
    # compute Monte-Carlo p-value
    mcb_resamples = sapply(1:n_resamples,function(i) score(x,resamples[,i]) - score(x_rc_resamples[,i],resamples[order(x),i]))
    rank_obs = tail(rank(c(mcb_resamples,mcb)),1)
    pval = 1 - (rank_obs - 1)/(n_resamples + 1)
  }

  text(x = 0,y = 1,paste0(c("BS ","MCB ","DSC ","UNC "),collapse = "\n"),adj = c(0,1))
  text(x = 0.2,y = 1,paste0(bquote(.(format(round(c(s,mcb,dsc,unc),digits = digits)),nsmall = digits)),collapse = "\n"),adj = c(0,1))
      
  abline(a = 0,b = 1,col = "grey",lty = 2)
  points(sort(x),x_rc,type = "l")

  if(inset_hist){
    a = par("usr")
    a = c(grconvertX(a[1:2], "user", "ndc"),grconvertY(a[3:4], "user", "ndc"))
    par.old = par(fig = c(0.3*a[1] + 0.7*a[2],0.05*a[1] + 0.95*a[2],0.9*a[3] + 0.1*a[4],0.65*a[3] + 0.35*a[4]),
                  pty = "m",mar = c(1,0,0,0),mgp = c(1,0.4,0),tcl = -0.25,new = TRUE)
    plot(hist(x,main = "",yaxt = "n",xlab = "",ylab = ""),add = TRUE)
    par(par.old)
  }
  
  return(list(decomp = c(mcb,dsc,unc),pval = if(resampling) pval else NULL))
}

# Mean and quantile reliability diagrams
reldiag = function(x,y,
                   type = "mean",      # for quantiles set type = list("quantile",alpha = 0.1) (replace 0.1 with appropriate level)
                   resampling = TRUE,n_resamples = 1000,replace = TRUE,region_level = 0.9,
                   digits = 3,inset_hist = TRUE,hist.breaks = 8,scatter_plot = TRUE,
                   lim = NULL,         # plotting region (used for xlim and ylim), e.g. lim = c(0,1)
                   main = "",xlab = NULL,ylab = NULL,adj_xlab = NA){
  if(type[[1]] == "mean"){
    pava = function(x,y) isoreg(x,y)$yf
    score = function(x,y) mean((x-y)^2) # canonical score
    marg = base::mean
    identif = function(x,y) x-y
    if(is.null(xlab)) xlab = expression(x == m[1](F))
    score_label = "MSE "
  }
  else{
    if(type[[1]] == "quantile"){
      require(isotone)
      alpha = type[[2]]
      pava = function(x,y){
        # In case of ties, isotone::gpava uses the conditional mean instead of quantile, try e.g.,
        # gpava(c(-1,-1,-1),c(-1,0,0),solver = weighted.median,ties = "secondary")
        # the following step replaces y values with the respective quantile in case of ties
        y = unlist(lapply(split(y,x),function(y) rep(quantile(y,alpha,type = 1),length(y))),use.names = FALSE)
        return(gpava(x,y,solver = weighted.fractile,p = alpha,ties = "secondary")$x)
      }
      score = function(x,y) mean(2*(as.numeric(x >= y) - alpha)*(x-y))
      marg = function(x) quantile(x,alpha,type = 1)
      identif = function(x,y) as.numeric(x > y) - alpha
      if(is.null(xlab)) xlab = bquote(x == q[.(alpha)](F))
      score_label = "QS "
    }
    else stop("type must be \"mean\" or list(\"quantile\",level)")
  }
  
  if(is.null(lim)){
    lim = range(x) + c(-1,1)*max(abs(range(x)))*0.2
  }
  
  if(is.null(ylab)) ylab = expression(hat(x)[rc])

  plot(NULL,xlim = lim,ylim = lim,main = main,xlab = "",ylab = ylab)
  mtext(xlab,side = 1,line = par()$mgp[1],cex = par()$cex,adj = adj_xlab)
  
  ord_x = order(x)
  x = x[ord_x]
  y = y[ord_x]
  
  x_rc = pava(x,y)
  
  # We encountered suboptimal solutions for quantiles in rare cases, e.g.,
  # gpava(c(3,3,2,1,1),1:5,solver = weighted.fractile,p = 0.75, ties = "secondary")
  # which led to slightly negative DSC components
  # applying gpava a second time seems to fix this:
  if(type[[1]] == "quantile"){
    x_rc2 = pava(x_rc,y)
    if(!all(x_rc == x_rc2)){
      warning("Encountered multiple gpava solutions...")
      x_rc = x_rc2
    }
  }
  
  res = y - x
  
  # score decomposition
  s = score(x,y)
  c_rc_ucond = optim(par = 0,fn = function(c) score(x+c,y),method = "Brent",lower = min(res),upper = max(res))$par
  s_rc_ucond = score(x + c_rc_ucond,y)
  s_rc = score(x_rc,y)
  s_mg = score(marg(y),y)
  
  mcb = s - s_rc
  umcb = s - s_rc_ucond
  cmcb = s_rc_ucond - s_rc
  dsc = s_mg - s_rc
  unc = s_mg
  
  # compute p-value for hypothesis of unconditional calibration (t-test)
  v = identif(x,y)
  t = sqrt(length(v)) * mean(v)/sd(v)
  pval_ucond = 1 - abs(pt(t,length(v)-1) - 0.5)*2

  if(resampling){
    low = floor(n_resamples * (1-region_level)/2)
    up = n_resamples - low

    resamples = sapply(1:n_resamples,function(i) x + sample(res,length(y),replace = replace)) 
    x_rc_resamples = apply(resamples, 2, function(y) pava(x,y))
    x_rc_resamples_sorted = apply(x_rc_resamples,1,sort) - marg(res)

    ran_x = range(x)
    polygon(c(ran_x[1],x,ran_x[2],rev(x),ran_x[1]),
            c(ran_x[1],x_rc_resamples_sorted[up,],ran_x[2],rev(x_rc_resamples_sorted[low,]),ran_x[1]),
            border = NA,col = "lightblue1")
    points(x,x_rc_resamples_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(x,x_rc_resamples_sorted[up,],type = "l",lty = 1,col = "lightblue2")
    box()
    
    # compute Monte-Carlo p-value (for hypothesis of conditional calibration)
    mcb_resamples = sapply(1:n_resamples,function(i) score(x,resamples[,i]) - score(x_rc_resamples[,i],resamples[,i]))
    mcb_bounds = sort(c(mcb,mcb_resamples))[c(low,up)]
    rank_obs = tail(rank(c(mcb_resamples,mcb)),1)
    pval = 1 - (rank_obs - 1)/(n_resamples + 1)
  }
  text_pos = legend("topleft",legend = parse(text = c(score_label,expression(MCB[u]),expression(MCB[c]),"DSC","UNC")),plot = FALSE)
  text(x = lim[1],y = text_pos$text$y,parse(text = c(score_label,expression(MCB[u]),expression(MCB[c]),"DSC","UNC")),adj = c(0,0.5))
  text(x = 0.8*lim[1] + 0.2*lim[2],y = text_pos$text$y,
       bquote(.(format(round(c(s,umcb,cmcb,dsc,unc),digits = digits)),nsmall = digits)),
       adj = c(0,0.5))

  abline(a = 0,b = 1,col = "grey",lty = 2)
  points(x,x_rc,type = "l")
  
  if(scatter_plot){
    points(x,y,pch = 20,col = adjustcolor("black",alpha = 0.25),cex = 0.5)
  }
  
  if(inset_hist){
    a = par("usr")
    a = c(grconvertX(a[1:2], "user", "ndc"),grconvertY(a[3:4], "user", "ndc"))
    par.old = par(fig = c(0.3*a[1] + 0.7*a[2],0.05*a[1] + 0.95*a[2],0.9*a[3] + 0.1*a[4],0.65*a[3] + 0.35*a[4]),
                  pty = "m",mar = c(1,0,0,0),mgp = c(1,0.4,0),tcl = -0.25,new = TRUE)
    plot(hist(x,breaks = hist.breaks,main = "",yaxt = "n",xlab = "",ylab = "",xlim = lim),add = TRUE)
    par(par.old)
  }
  
  return(list(x = x,y = y,res = res,x_rc = x_rc,
                 decomp = c(umcb,cmcb,mcb,dsc,unc),
                 pval_u = pval_ucond,pval_c = if(resampling) pval else NA))
}

#############################################################
# Simulation setup

# Simulate perfect, unfocused and lopsided forecast
setup_normal = function(n,tau0 = 1.5,eta0 = 0.7,set_seed = TRUE,seed = 100){
  if(set_seed) set.seed(seed)
  mu = rnorm(n)
  tau = sample(tau0*c(-1,1),n,replace = TRUE)
  eta = sample(eta0*c(-1,1),n,replace = TRUE)
  y = rnorm(n) + mu
  
  F_perf = function(x) pnorm(x,mu)
  F_unf = function(x) 0.5*(pnorm(x,mu) + pnorm(x,mu+tau))
  F_lop = function(x) ifelse(x <= mu,
                             (1-eta)*pnorm(x,mu),
                             (1+eta)*pnorm(x,mu)-eta)
  
  m_perf = cbind(mu, mu^2 + 1, mu^3 + 3*mu)
  m_unf = cbind(mu + 0.5*tau,
                mu^2 + tau*mu + 0.5*tau^2 + 1,
                mu^3 + 1.5*tau*mu^2 + 3*(0.5*tau^2 + 1)*mu + 0.5*(tau^3 + 3*tau))
  phi0 = dnorm(0)
  m_lop = cbind(mu + 2*eta*phi0,
                mu^2 + 1 + 4*eta*phi0*mu,
                mu^3 + 3*mu + 2*eta*phi0*(3*mu^2 + 2))
  
  sample_perf = function(cases = 1:n) rnorm(length(cases)) + mu[cases]
  sample_unf = function(cases = 1:n){
    mix = sample(c(0,1),length(cases),replace = TRUE)
    return(rnorm(length(cases)) + mu[cases] + mix*tau[cases])
  }
  sample_lop = function(cases = 1:n){
    mix = sample(c(-1,1),length(cases),replace = TRUE,prob = c(1-eta0,1+eta0))*sign(eta[cases])
    return(mu[cases] + mix*abs(rnorm(length(cases))))
  }
  
  quant_perf = function(alpha) qnorm(alpha,mean = mu)
  quant_unf = function(alpha) sapply(1:n, function(i) uniroot(function(x) 0.5*(pnorm(x,mu[i]) + pnorm(x,mu[i]+tau[i])) - alpha, 
                                                              c(qnorm(alpha,mu[i]-tau0), qnorm(alpha,mu[i]+tau0)))$root)
  quant_lop = function(alpha) ifelse(alpha <= (1-eta)/2, qnorm(pmin(1,alpha/(1-eta)),mu), qnorm(pmax(0,(alpha+eta)/(1+eta)),mu))
  
  return(list(mu = mu, tau = tau, eta = eta, y = y,
              perf = list(F = F_perf, quant = quant_perf, m = m_perf, sample = sample_perf),
              unf = list(F = F_unf, quant = quant_unf, m = m_unf, sample = sample_unf),
              lop = list(F = F_lop, quant = quant_lop, m = m_lop, sample = sample_lop)))
}

#############################################################
# Data example: BoE CPI inflation forecasts

# Distribution and sampling functions for Two-Piece-Normal distributions (Formulas from Julio (2006))
ptpnorm = function(x,mode,mean,unc){
  xi = mean - mode
  bet = pi/(2*unc^2)*xi^2
  gam = ifelse(xi == 0,0,
               -sign(xi)*sqrt(1 - ((sqrt(1 + 2*bet) - 1)/bet)^2))
  sig1 = unc/sqrt(1-gam)
  sig2 = unc/sqrt(1+gam)
  C = sqrt(2/pi)/(sig1 + sig2)
  ifelse(mode != mean,
         ifelse(x <= mode,
                C * sqrt(2*pi) * sig1 * pnorm(x,mode,sig1),
                1 - C * sqrt(2*pi) * sig2 * (1 - pnorm(x,mode,sig2))),
         pnorm(x,mode,unc))
}

qtpnorm = function(p,mode,mean,unc){
  xi = mean - mode
  bet = pi/(2*unc^2)*xi^2
  gam = ifelse(xi == 0,0,
               -sign(xi)*sqrt(1 - ((sqrt(1 + 2*bet) - 1)/bet)^2))
  sig1 = unc/sqrt(1-gam)
  sig2 = unc/sqrt(1+gam)
  C = sqrt(2/pi)/(sig1 + sig2)
  
  ifelse(p <= ptpnorm(mode,mode,mean,unc),
         mode + sig1*qnorm(p/(C*sqrt(2*pi)*sig1)),
         mode + sig2*qnorm((p+C*sqrt(2*pi)*sig2 - 1)/(C*sqrt(2*pi)*sig2)))
}

rtpnorm = function(mode,mean,unc){
  xi = mean - mode
  bet = pi/(2*unc^2)*xi^2
  gam = ifelse(xi == 0,0,
               -sign(xi)*sqrt(1 - ((sqrt(1 + 2*bet) - 1)/bet)^2))
  sig1 = unc/sqrt(1-gam)
  sig2 = unc/sqrt(1+gam)
  C = sqrt(2/pi)/(sig1 + sig2)
  
  n = length(mode)
  u = runif(n)
  
  p_mode = ptpnorm(mode,mode,mean,unc)
  
  ifelse(p_mode <= u,
         mode - abs(rnorm(n,sd = sig1)),
         mode + abs(rnorm(n,sd = sig2)))
}

process_data_CPI = function(raw_fcasts, raw_obs){
  quarters = raw_fcasts[1,-(1:4)]
  # Note: First two quarters of 2004 have less 'market' forecast values. Manually account for this at beginning of row indexes
  n = 68
  n_lags = 8
  rstart = 18
  cstart = 5

  # extract forecasts based on market interest
  rper10 = c(1,3,4)
  rows = c(7,9,5,15,17,13, rstart - 1 + rep(10*0:(n-3),each = 3) + rper10)
  
  # arrange forecasts in a data frame
  fcast_data = data.frame(Q = rep(unlist(quarters[1:n]),each = 3),par = raw_fcasts[rows,4])
  
  for(lag in 0:n_lags){
    cols = cstart - 1 + rep(1:n,each = 3) + lag
    assign(paste0("lag",lag),raw_fcasts[cbind(rows,cols)])
    fcast_data = data.frame(fcast_data,as.numeric(get(paste0("lag",lag))))
    names(fcast_data)[3 + lag] = paste0("lag",lag)
  }
  
  ind_first_obs = which(raw_obs[,1] == "2004 Q1")
  obs = data.frame(Q = unlist(quarters[1:n]),y = as.numeric(raw_obs[cbind(ind_first_obs + 0:(n-1),2)]))
  
  data_at_lag = function(lag = 1,max_lag = 6){
    max_obs = 65    # We use the first 65 observations, because the last 3 are missing nowcasts
    ind_obs = (max_lag+1):max_obs # We use those observations for which forecasts at horizons up to six quarters are available

    mode = fcast_data[paste0("lag",lag)][fcast_data$par == "Market Mode",][ind_obs - lag]
    mean = fcast_data[paste0("lag",lag)][fcast_data$par == "Market Mean",][ind_obs - lag]
    unc = fcast_data[paste0("lag",lag)][fcast_data$par == "Uncertainty",][ind_obs - lag]

    y = obs$y[ind_obs]
    F_fcast = function(x) ptpnorm(x,mode,mean,unc)
    m_fcast = cbind(mean,NA,NA)
    q_fcast = function(p) mapply(function(mode,mean,unc) qtpnorm(p,mode,mean,unc),mode,mean,unc)
    sample_fcast = function(cases = 1:length(y)) rtpnorm(mode[cases],mean[cases],unc[cases])
    
    return(list(fcast = list(F = F_fcast,quant = q_fcast,m = m_fcast,sample = sample_fcast),y = y,mode = mode, mean = mean,unc = unc))
  }
  return(data_at_lag)
}



