# The following file need to be in the working directory
# Warning: The 'load_data.R' script downloads some files to the working directory
# source("load_data.R")
# source("functions.R")

options(warn = -1)

mai_def = c(0.3,0.3,0.1,0.1)
mgp_def = c(-1.5,0.8,0)
lwd_def = 1.5
cex.main_def = 1.5
cex_def = 0.67

# Warning: A new subdirectory 'figs' is created and all figures will be saved there as png
# Create new subdirectory for figures
figs_path = "./figs/"
if(!dir.exists(figs_path)) dir.create(figs_path)

# Figure 2
plot_fig2 = function(file,width = 5.25,height = 2.7){
  cdf_vals = matrix(c(0,0,5,10,12,20,20,
                      0,0,5,10,18,20,20,
                      0,0,8,10,12,20,20,
                      0,0,2,10,18,20,20),ncol = 4)/20 + rep(0.01 * c(0.5,-0.5,1.5,-1.5),each = 7)
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 200)
  par(pty = "m",mai = mai_def,mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  matplot(c(-1,0,1,2,3,4,5),cdf_vals,type = "l",col = c("red","blue")[c(1,2,1,2)],lty = c(1,1,2,2),xlim = c(0,4),
          xlab = expression(y),ylab = expression(F(y)))
  legend("bottomright",
         c(expression(F[1]),expression(F[2]),expression(Q(Y<=y ~"|"~ F == F[1])),expression(Q(Y<=y ~"|"~ F == F[2]))),
         col = c("red","blue")[c(1,2,1,2)],lty = c(1,1,2,2))
  dev.off()
}

plot_fig2("fig2.png")

# Figure 3
plot_fig3 = function(plot_fct,fcast_fct,param,file,width,height,
                     mai = mai_def,main = "",main2 = "", ...){
  colors = c(brewer.pal(n = 8,name = "YlOrRd")[5:8],
             "black",
             brewer.pal(n = 8,name = "YlGnBu")[8:5])
  png(file = paste0(figs_path,file),width = width, height = height,units = "in",res = 200)
  par(pty = "s",mai = mai,mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  plot_fct(fcast_fct,param,main = "",colors = colors, ...)
  title(main = main,line = 1)
  title(ylab = main2,font.lab = par()$font.main,cex.lab = par()$cex.main,line = 2.8)
  dev.off()
}

plot_fig3(plot_thresh,cond_nep_unf,-2:2,"fig3a.png",width = 3.15,height = 3.15,
          mai = mai_def + c(0,0.45,0.45,0),main = "Unfocused",main2 = "Threshold")
plot_fig3(plot_thresh,cond_nep_lop,-2:2,"fig3b.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Lopsided",main2 = "")
plot_fig3(plot_thresh,cond_nep_pwUnif,seq(0,3,0.75),"fig3c.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Piecewise Uniform",main2 = "")
plot_fig3(plot_quant,cond_quant_unf,seq(0.1,0.9,0.1),"fig3d.png",width = 3.15,height = 2.7,
          mai = mai_def + c(0,0.45,0,0),main = "",main2 = "Quantile")
plot_fig3(plot_quant,cond_quant_lop,seq(0.1,0.9,0.1),"fig3e.png",width = 2.7,height = 2.7)
plot_fig3(plot_quant,cond_quant_pwUnif,seq(0.1,0.9,0.1),"fig3f.png",width = 2.7,height = 2.7,lim = c(-1,4))
plot_fig3(plot_moms,cond_moms_unf,NULL,"fig3g.png",width = 3.15,height = 2.7,
          mai = mai_def + c(0,0.45,0,0),main = "",main2 = "Moment")
plot_fig3(plot_moms,cond_moms_lop,NULL,"fig3h.png",width = 2.7,height = 2.7,mai = mai_def)
plot_fig3(plot_moms,cond_moms_pwUnif,NULL,"fig3i.png",width = 2.7,height = 2.7,
          mai = mai_def,lim = c(-2,5))

# Figure S1
# uses plot_fig3 from Figure 3
plot_fig3(plot_moms,cond_moms_unf,NULL,"figS1a.png",width = 3.15,height = 3.15,
          mai = mai_def + c(0,0.45,0.45,0),main = "Unfocused",main2 = "Moment", lim = c(-64,64),root_transform = FALSE)
plot_fig3(plot_moms,cond_moms_lop,NULL,"figS1b.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Lopsided", lim = c(-64,64),root_transform = FALSE)
plot_fig3(plot_moms,cond_moms_pwUnif,NULL,"figS1c.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Piecewise Uniform", lim = c(-2,2),root_transform = FALSE)


# Figure 4
# Expected score of unfocused mean
plot_fig4 = function(file,width = 5.25,height = 2.7){
  Psi = function(x,c = 1) dnorm(x + c)/(dnorm(x) + dnorm(x + c))
  I = function(c) integrate(function(x) Psi(x,c)^2*dnorm(x),lower = -20,upper = 20)$value
  mcb = function(c) c^2*(0.25 - I(c))
  dsc = function(c) 1 - I(c)*c^2
  unc = 2
  
  c = seq(0,5,0.1)
  
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 200)
  par(pty = "m",mai = mai_def + c(0,0,0.45,0),mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  
  plot(c,rep(unc,length(c)),ylim = c(0,3),col = "blue",type = "l",xlab = expression(eta[0]),ylab = "",xaxs = "i")
  points(c,sapply(c,mcb),type = "l",col = "red")
  points(c,sapply(c,dsc),type = "l",col = "green")     
  points(c, sapply(c,function(c) mcb(c) - dsc(c) + unc),type = "l")
  legend("topright",legend = c("MSE","MCB","DSC","UNC"),lty = 1,col = c("black","red","green","blue"),bg = "white")
  
  labels = seq(0,1,0.5)
  at = labels*sqrt(8/pi)
  axis(3, at = at,labels = labels)
  mtext(expression(delta[0]), 3 ,line = -1.5,at = at[2],cex = cex_def)
  abline(v = at[3],col = "grey",lty = 2)
  dev.off()
}

plot_fig4("fig4.png")


# Figure 5
plot_fig5 = function(plot_reldiag,fcast,y,file,width,height,
                     mai = mai_def,main = "",main2 = "",...){
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 300)
  par(pty = "s",mai = mai,mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  plot_reldiag(fcast,y, ...)
  title(main = main,line = 1)
  title(ylab = main2,font.lab = par()$font.main,cex.lab = par()$cex.main,line = 2.8)
  dev.off()
}

sim = setup_normal(400,1.5,0.7)

plot_fig5(threshreldiag,sim$perf,sim$y,"fig5a.png",width = 3.15,height = 3.15,
          mai = mai_def + c(0,0.45,0.45,0),main = "Perfect",main2 = "Threshold",t = 1)
plot_fig5(threshreldiag,sim$unf,sim$y,"fig5b.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Unfocused",main2 = "",t = 1)
plot_fig5(threshreldiag,sim$lop,sim$y,"fig5c.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Lopsided",main2 = "",t = 1)

lim_mean = c(-4.5,4.5)
plot_fig5(reldiag,sim$perf$m[,1],sim$y,"fig5d.png",width = 3.15,height = 2.7,
          mai = mai_def + c(0,0.45,0,0),main = "",main2 = "Mean",type = "mean",scatter_plot = TRUE,inset_hist = FALSE, lim = lim_mean)
plot_fig5(reldiag,sim$unf$m[,1],sim$y,"fig5e.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "",type = "mean",scatter_plot = TRUE,inset_hist = FALSE, lim = lim_mean)
plot_fig5(reldiag,sim$lop$m[,1],sim$y,"fig5f.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "",type = "mean",scatter_plot = TRUE,inset_hist = FALSE, lim = lim_mean)

lim_quant = c(-6,4)
plot_fig5(reldiag,sim$perf$quant(0.1),sim$y,"fig5g.png",width = 3.15,height = 2.7,
          mai = mai_def + c(0,0.45,0,0),main = "",main2 = "Quantile",type = list("quantile",alpha = 0.1),scatter_plot = FALSE, lim = lim_quant)
plot_fig5(reldiag,sim$unf$quant(0.1),sim$y,"fig5h.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "",type = list("quantile",alpha = 0.1),scatter_plot = FALSE, lim = lim_quant)
plot_fig5(reldiag,sim$lop$quant(0.1),sim$y,"fig5i.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "",type = list("quantile",alpha = 0.1),scatter_plot = FALSE, lim = lim_quant)

# Figure 7
# uses plot_fig5 from Figure 5!
plot_fig5(reldiag,obs_pred$pred_null,obs_pred$obs,"fig7a.png",width = 4,height = 4.45,
          mai = mai_def + c(0,0,0.45,0),main = "Null Model",main2 = "",
          type = "mean",lim = c(-0.35,3),scatter_plot = TRUE,inset_hist = FALSE,adj_xlab = 0.9,)
plot_fig5(reldiag,obs_pred$pred_opt,obs_pred$obs,"fig7b.png",width = 4,height = 4.45,
          mai = mai_def + c(0,0,0.45,0),main = "Ridge Regression",main2 = "",
          type = "mean",lim = c(-0.35,3),scatter_plot = TRUE,inset_hist = FALSE,adj_xlab = 0.9)

# Figure 8 and S3 - S8
# uses plot_fig5 for reliability diagrams and plot_fig8 for ACFs
plot_fig8 = function(z,file,width,height,mai = mai_def,main = "",main2 = "",...){
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 200)
  par(pty = "s",mai = mai,mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  acf(z,main = "",xlab = "lag", ...)
  title(main = main,line = 1)
  title(ylab = main2,font.lab = par()$font.main,cex.lab = par()$cex.main,line = 2.8)
  dev.off()
}

data_at_lag = process_data_CPI(raw_fcasts,raw_obs)

# lim_mean = c(-1,6)
lim_quant = c(0,6)
fig_labels = c("S3","8","S4","S5","S6","S7","S8")

for(lag in 0:6){
  data_eval = data_at_lag(lag)
  z = data_eval$fcast$F(data_eval$y)
  z2 = (z - 0.5)^2
  
  plot_fig5(PITreldiag,data_eval$fcast,data_eval$y,paste0("fig",fig_labels[lag+1],"a.png"),width = 2.7,height = 3.15,
            mai = mai_def + c(0,0,0.45,0),main = "(a) PIT",main2 = "")
  plot_fig8(z,paste0("fig",fig_labels[lag+1],"b.png"),width = 2.7,height = 3.15,
                mai = mai_def + c(0,0,0.45,0),main = "(b) ACF of PIT",main2 = "",ci = 0.9,ylab = "")
  plot_fig8(z2,paste0("fig",fig_labels[lag+1],"c.png"),width = 2.7,height = 3.15,
                mai = mai_def + c(0,0,0.45,0),main = expression(bold("(c) ACF of" ~ group("(",PIT - 1/2,")")^2)),main2 = "",ci = 0.9,ylab = "")
  plot_fig5(margreldiag,data_eval$fcast,data_eval$y,paste0("fig",fig_labels[lag+1],"d.png"),width = 2.7,height = 3.15,
            mai = mai_def + c(0,0,0.45,0),main = "(d) Marginal",main2 = "")
  plot_fig5(threshreldiag,data_eval$fcast,data_eval$y,paste0("fig",fig_labels[lag+1],"e.png"),width = 2.7,height = 3.15,
            mai = mai_def + c(0,0,0.45,0),main = "(e) Threshold",main2 = "",t = 2)
  plot_fig5(reldiag,data_eval$fcast$quant(0.75),data_eval$y,paste0("fig",fig_labels[lag+1],"f.png"),width = 2.7,height = 3.15,
            mai = mai_def + c(0,0,0.45,0),main = "(f) Quantile",main2 = "",type = list("quantile",0.75),scatter_plot = FALSE,lim = lim_quant)
}

# Figure S2
# uses plot_fig5 from Figure 5!
plot_fig5(PITreldiag,sim$perf,sim$y,"figS2a.png",width = 3.15,height = 3.15,
          mai = mai_def + c(0,0.45,0.45,0),main = "Perfect",main2 = "PIT")
plot_fig5(PITreldiag,sim$unf,sim$y,"figS2b.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Unfocused",main2 = "")
plot_fig5(PITreldiag,sim$lop,sim$y,"figS2c.png",width = 2.7,height = 3.15,
          mai = mai_def + c(0,0,0.45,0),main = "Lopsided",main2 = "")
plot_fig5(margreldiag,sim$perf,sim$y,"figS2d.png",width = 3.15,height = 2.7,
          mai = mai_def + c(0,0.45,0,0),main = "",main2 = "Marginal")
plot_fig5(margreldiag,sim$unf,sim$y,"figS2e.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "")
plot_fig5(margreldiag,sim$lop,sim$y,"figS2f.png",width = 2.7,height = 2.7,
          mai = mai_def,main = "",main2 = "Marginal")


# Figure 6
plot_fig6 = function(file,width = 5,height = 2.7){
  x = c(1,2,4,6,8,10,11,12,14)
  y = c(4,5,6,9,10,11,13,8,15)
  
  olsreg = lm(y~x)
  
  require(quantreg)
  medreg = rq(y~x)
  
  m1_rc = isoreg(y~x)
  rep(0.01 * c(0.5,-0.5,1.5,-1.5),each = 7)

  # Isotonic Median Regression
  require(isotone)
  med_rc = gpava(x,y,solver = weighted.median)$x
  med_rc_up = med_rc
  med_rc_up[7:8] = 13
  
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 200)
  
  # plotting adjustments:
  m1_rc_y = m1_rc$yf + rep(-0.05,length(x))
  med_rc = med_rc + rep(0.05,length(x))
  med_rc_up = med_rc_up + rep(0.05,length(x))
  
  par(pty = "s",mai = mai_def + c(0,0,0,2.2),mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)

  plot(x,y,xlim = c(0,15),ylim = c(2,17),xlab = "x",ylab = "y")
  
  polygon(c(x[6:9],x[8:7]),c(med_rc_up[6:9],med_rc[8:7]),col = adjustcolor("green",alpha = 0.15),border = NA)
  
  abline(olsreg$coefficients,col = "blue")
  abline(medreg$coefficients,col = "green")

  points(x, med_rc,type = "l",col = "green",lty = 2)
  points(x, med_rc_up,type = "l",col = "green",lty = 2)
  points(x,m1_rc_y,type = "l",col = "blue",lty = 2)
  
  par(xpd = TRUE)
  
  coeff_det = function(x,type = "mean"){
    mse = function(x,y) sum((x-y)^2)/length(y)
    mae = function(x,y) sum(abs(x-y))/length(y)
    if(type == "mean"){
      x_rc = m1_rc$yf
      score = mse
    }
    else if(type == "median"){
      x_rc = med_rc
      score = mae
    }
    else stop("type must be mean or median!")
    mean_score = score(x,y)
    mean_score_rc = score(x_rc,y)
    mean_score_mg = score(mean(y),y)
    
    return((mean_score_mg - mean_score)/mean_score_mg)
  }

  legend("bottomright", inset = c(-0.8,0),
              legend = c("Linear mean regression",
                         as.expression(bquote(list(R^2 == .(format(round(coeff_det(olsreg$fitted.values),3),nsmall = 3)),
                                     R^1 == .(format(round(coeff_det(olsreg$fitted.values,"median"),3),nsmall = 3))))),
                         "Linear median regression",
                         as.expression(bquote(list(R^2 == .(format(round(coeff_det(medreg$fitted.values),3),nsmall = 3)),
                                     R^1 == .(format(round(coeff_det(medreg$fitted.values,"median"),3),nsmall = 3))))),
                         "Isotonic mean regression","Isotonic median regression"),
              col = c("blue","green")[c(1,NA,2,NA,1,2)],lty = c(1,NA,1,NA,2,2),bg = "white")
 dev.off()
}

plot_fig6("fig6.png")

# Figure 9
plot_fig9 = function(mdu,label_score = "MSE",file,width,height,
                     mai = mai_def,main = "",show_MCBu = FALSE,...){
  png(file = paste0(figs_path,file),width = width,height = height,units = "in",res = 200)
  layout(matrix(c(1,2),nrow = 2),heights = c(2, 1))
  par(pty = "m",mai = mai,mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)

  cols = c("black","red","green","blue","black")
  lags = 0:6
  
  layout(matrix(c(1,2),nrow = 2),heights = c(2, 1))
  par(pty = "m",mai = c(0,mai[2:4]),mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  
  if(show_MCBu) ylim = c(0,1.3*max(mdu[,1:4]))
  else ylim = c(0,1.3*max(mdu[,1:4]))
  
  matplot(lags,mdu[,1:4],type = "b",
          pch = 4,col = cols[1:4],lty = rep(1,4),ylim = ylim, xaxt = "n", xlab = "",ylab = "")
  if(show_MCBu) points(lags,mdu[,6],type = "b",pch = 4,col = "red",lty = 2)
  abline(h = 0,col = "grey",lty = 2)
  if(show_MCBu) legend("topleft",legend = parse(text = c(label_score,"MCB","DSC",expression(MCB[u]),"UNC")),lty = c(1,1,1,2,1),col = c(cols[c(1,2,3,2,4)]),ncol = 3,bg = 'white')
  else legend("topleft",legend = c(label_score,"MCB","DSC","UNC"),lty = 1,col = cols[1:4],ncol = 2,bg = 'white')
  title(main = main,line = 1)

  par(pty = "m",mai = c(mai[1:2],0.05,mai[4]),mgp = mgp_def,lwd = lwd_def,cex = cex_def,cex.main = cex.main_def)
  plot(lags,mdu[,5],type = "b",ylim = c(-0.5,1),
          pch = 4,col = cols[5],lty = 2,xlab = "lead time",ylab = "")
  abline(h = 0,col = "grey",lty = 2)
  legend("topright",legend = expression(R^x),lty = 2,col = cols[5])

  dev.off()
}

data_at_lag = process_data_CPI(raw_fcasts,raw_obs)

decomp_quant = function(alpha){
  mdu = sapply(0:6,function(lag){
    data_eval = data_at_lag(lag)
    return(reldiag(data_eval$fcast$quant(alpha),data_eval$y,type = list("quantile",alpha),resampling = FALSE)$decomp)
  })
  return(cbind(mdu[3,] - mdu[4,] + mdu[5,],t(mdu[3:5,]),(mdu[4,]-mdu[3,])/mdu[5,],mdu[1,]))
}
decomp_quant = decomp_quant(0.75)

decomp_thresh = function(t = 2){
  mdu = sapply(0:6,function(lag){
    data_eval = data_at_lag(lag)
    return(threshreldiag(data_eval$fcast,data_eval$y,t,resampling = FALSE)$decomp)
  })
  return(cbind(mdu[1,] - mdu[2,] + mdu[3,],t(mdu),(mdu[2,]-mdu[1,])/mdu[3,]))
}
decomp_thresh = decomp_thresh()

plot_fig9(decomp_thresh,"BS","fig9a.png",width = 4,height = 3.45,
          mai = mai_def + c(0,0,0.45,0),mgp = c(-1.2,0.8,0),main = "Threshold")
plot_fig9(decomp_quant,"QS","fig9b.png",width = 4,height = 3.45,
          mai = mai_def + c(0,0,0.45,0),mgp = c(-1.2,0.8,0),main = "Quantile",show_MCBu = TRUE)


