corrMat = 0.6^as.matrix(dist(1:6))
V1 = diag(sqrt(xs[1,-1]))
V2 = diag(sqrt(xs[2,-1]))
Sigma1 = V1%*%corrMat%*%V1
Sigma2 = V2%*%corrMat%*%V2



runSimComp = function(alpha, n.iter, nx, ny, corr = 0.8, sigma1, sigma2){

  dev.y.1 = 2.21*alpha
  dev.y.2 = 5.38*alpha

  xm = -t(cbind(c(0, 0.738, 1.313, 3.109, 4.525, 5.864, 7.338),
                c(0, 0.668, 3.975, 5.931, 8.259, 11.989, 13.976)))

  xs = t(cbind(c(0, 4.79, 5.43, 6.54, 7.37, 8.15, 9.11),
               c(0, 10.27, 12.85, 14.95, 15.35, 16.87, 18.19)))

  ym = xm + t(cbind(c(0,(rep(dev.y.1/6,6))), c(0,(rep(dev.y.2/6,6)))))

  lmm.sim = numeric(n.iter)
  lrst.sim = numeric(n.iter)
  nlpd.sim = numeric(n.iter)

  lrst.discrete = numeric(n.iter)
  nlpd.discrete = numeric(n.iter)


  for(iter in 1:n.iter){

    X_placebo = array(NA, dim = c(nx, 2, 7))
    X_dose_1 = array(NA, dim = c(ny, 2, 7))

    for(i in 1:nx){
      X_placebo[i,,] = gen.fun.dep(xm, xs, corr, sigma1, sigma2)
    }

    for(i in 1:ny){
      X_dose_1[i,,] = gen.fun.dep(ym, xs, corr, sigma1, sigma2)
    }


    X = X_placebo
    Y = X_dose_1
    K = dim(X)[2]
    T = dim(X)[3]-1

    X_c=array(0,dim=c(nx,K,T))
    for (i in 1:nx){for (j in 1:T) {X_c[i,,j]=X[i,,(j+1)]-X[i,,1]}}
    Y_c=array(0,dim=c(ny,K,T))
    for (i in 1:ny){for (j in 1:T) {Y_c[i,,j]=Y[i,,(j+1)]-Y[i,,1]}}


    lmm.sim[iter] = min(lmmFun4(X_c, Y_c))
    lrst.sim[iter] = uniUstat(X_c, Y_c)
    nlpd.sim[iter] = min(npld(X_c, Y_c))

    X_c_d = discretize(X_c)
    Y_c_d = discretize(Y_c)

    lrst.discrete[iter] = uniUstat(X_c_d, Y_c_d)
    nlpd.discrete[iter] = min(npld(X_c_d, Y_c_d))



    print(iter)
  }

  # print(sel.uni)
  lmS = length(which(lmm.sim<0.025))/length(lmm.sim)
  uss = length(which(lrst.sim<0.05))/length(lrst.sim)
  lds = length(which(nlpd.sim<0.025))/length(nlpd.sim)

  usd = length(which(lrst.discrete<0.05))/length(lrst.discrete)
  ldd = length(which(nlpd.discrete<0.025))/length(nlpd.discrete)

  return(list("LMM" = lmS, "LRST" = uss, "NPLD" = lds,
              "LRST.disc" = usd, "NPLD.disc" = ldd))
}

nx = 311   # Placebo Size
ny = 448   # Treatment Size
r = runSimComp(effect_Size, 1000, nx, ny, 0.8)
# Do the above for different effect sizes and different correlation



#####################   Outputs   #######################################################




alpha = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2)

power.LMM.0 = c(0.067, 0.121, 0.164, 0.229, 0.358, 0.544, 0.688, 0.855, 0.991)
power.LRST.0 = c(0.138, 0.271, 0.347, 0.452, 0.624, 0.791, 0.895, 0.976, 0.998)
power.NPLD.0 = c(0.064, 0.119, 0.154, 0.214, 0.352, 0.502, 0.665, 0.843, 0.983)
power.LRST.disc.0 = c(0.13, 0.261, 0.345, 0.433, 0.601, 0.77, 0.88, 0.968, 0.998)
power.NPLD.disc.0 = c(0.057, 0.122, 0.153, 0.203, 0.342, 0.463, 0.636, 0.817, 0.98)


power.LMM.2 = c(0.074, 0.121, 0.159, 0.195, 0.326, 0.524, 0.66, 0.868, 0.983)
power.LRST.2 = c(0.123, 0.235, 0.335, 0.384, 0.564, 0.753, 0.854, 0.962, 0.999)
power.NPLD.2 = c(0.065, 0.115, 0.159, 0.184, 0.31, 0.499, 0.639, 0.823, 0.981)
power.LRST.disc.2 = c(0.122, 0.231, 0.316, 0.366, 0.557, 0.74, 0.857, 0.957, 0.998)
power.NPLD.disc.2 = c(0.065, 0.106, 0.145, 0.171, 0.303, 0.484, 0.614, 0.812, 0.977)



power.LMM.5 = c(0.075, 0.133, 0.154, 0.205, 0.314, 0.468, 0.623, 0.831, 0.972)
power.LRST.5 = c(0.126, 0.221, 0.294, 0.358, 0.505, 0.686, 0.802, 0.935, 0.994)
power.NPLD.5 = c(0.077, 0.129, 0.146, 0.201, 0.304, 0.461, 0.604, 0.813, 0.964)
power.LRST.disc.5 = c(0.121, 0.218, 0.284, 0.35, 0.491, 0.663, 0.795, 0.935, 0.989)
power.NPLD.disc.5 = c(0.063, 0.132, 0.146, 0.185, 0.276, 0.447, 0.564, 0.796, 0.957)


power.LMM.8 = c(0.057, 0.112, 0.152, 0.183, 0.293, 0.458, 0.55, 0.782, 0.959)
power.LRST.8 = c(0.112, 0.219, 0.261, 0.318, 0.478, 0.626, 0.72, 0.887, 0.984)
power.NPLD.8 = c(0.065, 0.106, 0.148, 0.182, 0.284, 0.443, 0.535, 0.752, 0.952)
power.LRST.disc.8 = c(0.108, 0.213, 0.255, 0.314, 0.466, 0.612, 0.71, 0.871, 0.98)
power.NPLD.disc.8 = c(0.068, 0.106, 0.138, 0.172, 0.263, 0.441, 0.507, 0.727, 0.944)



df0 = data.frame(Deviation = c(alpha, alpha, alpha), Power = c(power.LMM.0, power.NPLD.0, power.LRST.0),
                 Method = c(rep("LMM", length(alpha)), rep("NPLD", length(alpha)), rep("LRST", length(alpha))))

p1 = ggplot(df0, aes(Deviation, Power, col = Method)) + geom_line(linewidth = 1)+ geom_point()+
  theme(axis.text.x = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.1, vjust = 0.2)) +
  scale_x_continuous("Deviation", labels = as.character(alpha), breaks = alpha)+
  # expand_limits(x=alpha)+
  theme(axis.text.y = element_text(colour = 'black', size = 12), axis.title.y = element_text(size = 12,
                                                                                             hjust = 0.5, vjust = 0.2))  + theme_classic()


df2 = data.frame(Deviation = c(alpha, alpha, alpha), Power = c(power.LMM.2, power.NPLD.2, power.LRST.2),
                 Method = c(rep("LMM", length(alpha)), rep("NPLD", length(alpha)), rep("LRST", length(alpha))))

p2 = ggplot(df2, aes(Deviation, Power, col = Method)) + geom_line(linewidth = 1)+ geom_point()+
  theme(axis.text.x = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.1, vjust = 0.2)) +
  scale_x_continuous("Deviation", labels = as.character(alpha), breaks = alpha)+
  # expand_limits(x=alpha)+
  theme(axis.text.y = element_text(colour = 'black', size = 12), axis.title.y = element_text(size = 12,
                                                                                             hjust = 0.5, vjust = 0.2))  + theme_classic()



df5 = data.frame(Deviation = c(alpha, alpha, alpha), Power = c(power.LMM.5, power.NPLD.5, power.LRST.5),
                 Method = c(rep("LMM", length(alpha)), rep("NPLD", length(alpha)), rep("LRST", length(alpha))))

p5 = ggplot(df5, aes(Deviation, Power, col = Method)) + geom_line(linewidth = 1)+ geom_point()+
  theme(axis.text.x = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.1, vjust = 0.2)) +
  scale_x_continuous("Deviation", labels = as.character(alpha), breaks = alpha)+
  # expand_limits(x=alpha)+
  theme(axis.text.y = element_text(colour = 'black', size = 12), axis.title.y = element_text(size = 12,
                                                                                             hjust = 0.5, vjust = 0.2))  + theme_classic()



df8 = data.frame(Deviation = c(alpha, alpha, alpha), Power = c(power.LMM.8, power.NPLD.8, power.LRST.8),
                 Method = c(rep("LMM", length(alpha)), rep("NPLD", length(alpha)), rep("LRST", length(alpha))))

p8 = ggplot(df8, aes(Deviation, Power, col = Method)) + geom_line(linewidth = 1)+ geom_point()+
  theme(axis.text.x = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.1, vjust = 0.2)) +
  scale_x_continuous("Deviation", labels = as.character(alpha), breaks = alpha)+
  # expand_limits(x=alpha)+
  theme(axis.text.y = element_text(colour = 'black', size = 12), axis.title.y = element_text(size = 12,
                                                                                             hjust = 0.5, vjust = 0.2))  + theme_classic()
par(mfrow= c(2,2))
gridExtra::grid.arrange(p1,p2,p5,p8, common.legend = TRUE, legend = "bottom")
library(ggpubr)
ggarrange(p1, p2, p5, p8, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


N = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000)
nxN = 2/5 * N
nyN = 3/5 * N
lmm.0 = numeric(13)
lrst.0 = numeric(13)
npld.0 = numeric(13)
lrst.disc = numeric(13)
npld.disc = numeric(13)
for(i in 3:6){
  nx = nxN[i]
  ny = nyN[i]
  r = runSimComp(0, 1000, nx, ny, 0.5, Sigma1, Sigma2)
  lmm.0[i] = r$LMM
  lrst.0[i] = r$LRST
  npld.0[i] = r$NPLD
  lrst.disc[i] = r$LRST.disc
  npld.disc[i] = r$NPLD.disc
  print(i)
}

lmm.0 = c(0.056, 0.047, 0.050, 0.047, 0.049, 0.049, 0.053, 0.048, 0.047, 0.041, 0.043, 0.048)
lrst.0 = c(0.058, 0.057, 0.047, 0.053, 0.050, 0.053, 0.052, 0.049, 0.053,0.059, 0.055, 0.053)
npld.0 = c(0.050, 0.043, 0.044, 0.050, 0.048, 0.046, 0.048, 0.049, 0.045 ,0.040, 0.037, 0.049)
lrst.disc.0 = c(0.055, 0.055, 0.050, 0.056, 0.054, 0.053, 0.054, 0.060, 0.051,0.054, 0.052, 0.051)
npld.disc.0 = c(0.053, 0.043, 0.045, 0.045, 0.051, 0.037, 0.051, 0.054, 0.049, 0.049, 0.033, 0.047)
