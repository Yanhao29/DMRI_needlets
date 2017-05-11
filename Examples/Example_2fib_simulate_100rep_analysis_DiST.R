rm(list=ls())
library(R.matlab)

path_load = '/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/'

sep = 45
lmax = 16
b = 1
ratio = 10
n_sample = 41

folder_name = paste0('2fib_sep',toString(sep),'_lmax',toString(lmax),'_b',toString(b),
                     '_ratio',toString(ratio),'_n',toString(n_sample),'_sig0.05/')

file_name = paste0('2fib_sep',toString(sep),'_lmax',toString(lmax),'_b',toString(b),
                     '_ratio',toString(ratio),'_n',toString(n_sample),'_sig0.05_rep')


source('/Users/hao/Dropbox/DMRI_code/DiST/dwi_basic.R')
source('/Users/hao/Dropbox/DMRI_code/DiST/dwi_fit.R')

opt = NULL ## set options for Raymonds algorithm
if (is.null(opt$maxJ)){
  opt$maxJ <- 4
}
if (is.null(opt$basis.depth)){
  opt$basis.depth <- 3
}
if (is.null(opt$basis.random)){
  opt$basis.random <- T
}
if (is.null(opt$def.alphas)){
  opt$def.alphas <- 2
}
if (is.null(opt$Ninner)){
  opt$Ninner <- 10
}
if (is.null(opt$Nmiddle.max)){
  opt$Nmiddle.max <- 50
}
if (is.null(opt$pre.thres)){
  opt$pre.thres <- 0.2
}

N_rep = 100

nfib_DiST = rep(0,N_rep)
angle1_DiST = rep(0,N_rep)
angle2_DiST = rep(0,N_rep)
sep_DiST = rep(0,N_rep)

theta1_DiST = rep(0,N_rep)
phi1_DiST = rep(0,N_rep)

theta2_DiST = rep(0,N_rep)
phi2_DiST = rep(0,N_rep)

fibVec_DiST = array(0,dim=c(N_rep,2,3))
fib_true_DiST = array(0,dim=c(N_rep,2,3))

fit_DiST = rep(list(),N_rep)

for(rep in 1:N_rep){
  ## load matlab simulation
  temp = readMat(paste0(path_load,folder_name,file_name,toString(rep),'.mat'))
  grad.mat_temp <- t(temp$pos.sampling)
  design_temp <- generate.design.noS0(grad.mat_temp)
  
  print(rep)
  
  out <- try(biased.mle.cs2.iso(dwis=c(temp$DWI), design=design_temp,
                                sigma=c(temp$sig), S0=1, b=c(temp$b)[1], maxJ=opt$maxJ,
                                basis.depth=opt$basis.depth,
                                basis.random=opt$basis.random,
                                def.alphas=opt$def.alphas, betas0=NULL,
                                Ninner=opt$Ninner, Nmiddle.max=opt$Nmiddle.max,
                                display=F, pre.thres=opt$pre.thres,
                                alphas.const=opt$alphas.const))
  nfib_DiST[rep] = out$n.fiber
  
  if(out$n.fiber == 2){
    fib_temp = out$vs
    theta.r = temp$theta.r - t(temp$theta0)
    phi.r = temp$phi.r - t(temp$phi0)
    
    fibVec_DiST[rep,,] = out$vs
    phi1_DiST[rep] = temp$phi.r[1]
    phi2_DiST[rep] = temp$phi.r[2]
    theta1_DiST[rep] = temp$theta.r[1]
    theta2_DiST[rep] = temp$theta.r[2]
    rotationM1 = matrix(c(cos(phi.r[1])*cos(theta.r[1]),-sin(phi.r[1]),
                          cos(phi.r[1])*sin(theta.r[1]),sin(phi.r[1])*cos(theta.r[1]),
                          cos(phi.r[1]),sin(phi.r[1])*sin(theta.r[1]),-sin(theta.r[1]),
                          0,cos(theta.r[1])),nrow=3)
    rotationM2 = matrix(c(cos(phi.r[2])*cos(theta.r[2]),-sin(phi.r[2]),
                          cos(phi.r[2])*sin(theta.r[2]),sin(phi.r[2])*cos(theta.r[2]),
                          cos(phi.r[2]),sin(phi.r[2])*sin(theta.r[2]),-sin(theta.r[2]),
                          0,cos(theta.r[2])),nrow=3)
    
    fib1_true = temp$fod1.s%*%rotationM1
    fib2_true = temp$fod2.s%*%rotationM2
    
    fib_true_DiST[rep,1,] = fib1_true
    fib_true_DiST[rep,2,] = fib2_true
    
    ang11 = min(acos(sum(fib_temp[1,]*fib1_true)/(sqrt(sum(fib_temp[1,]*fib_temp[1,]))*sqrt(sum(fib1_true*fib1_true)))),
               abs(pi-acos(sum(fib_temp[1,]*fib1_true)/(sqrt(sum(fib_temp[1,]*fib_temp[1,]))*sqrt(sum(fib1_true*fib1_true))))))*180/pi
    ang12 = min(acos(sum(fib_temp[2,]*fib1_true)/(sqrt(sum(fib_temp[2,]*fib_temp[2,]))*sqrt(sum(fib1_true*fib1_true)))),
                abs(pi-acos(sum(fib_temp[2,]*fib1_true)/(sqrt(sum(fib_temp[2,]*fib_temp[2,]))*sqrt(sum(fib1_true*fib1_true))))))*180/pi
    
    angle1_DiST[rep] = min(ang11,ang12)
    idx_dele = which.min(c(ang11,ang12))
    angle2_DiST[rep] = min(acos(sum(fib_temp[-idx_dele,]*fib2_true)/(sqrt(sum(fib_temp[-idx_dele,]*fib_temp[-idx_dele,]))*sqrt(sum(fib2_true*fib2_true)))),
                abs(pi-acos(sum(fib_temp[-idx_dele,]*fib2_true)/(sqrt(sum(fib_temp[-idx_dele,]*fib_temp[-idx_dele,]))*sqrt(sum(fib2_true*fib2_true))))))*180/pi
    sep_DiST[rep] = min(acos(sum(fib_temp[1,]*fib_temp[2,])/(sqrt(sum(fib_temp[1,]*fib_temp[1,]))*sqrt(sum(fib_temp[2,]*fib_temp[2,])))),
                        abs(pi-acos(sum(fib_temp[1,]*fib_temp[2,])/(sqrt(sum(fib_temp[1,]*fib_temp[1,]))*sqrt(sum(fib_temp[2,]*fib_temp[2,]))))))*180/pi
  } else {
    fib_temp = out$vs
    theta.r = temp$theta.r - t(temp$theta0)
    phi.r = temp$phi.r - t(temp$phi0)
    rotationM1 = matrix(c(cos(phi.r[1])*cos(theta.r[1]),-sin(phi.r[1]),
                          cos(phi.r[1])*sin(theta.r[1]),sin(phi.r[1])*cos(theta.r[1]),
                          cos(phi.r[1]),sin(phi.r[1])*sin(theta.r[1]),-sin(theta.r[1]),
                          0,cos(theta.r[1])),nrow=3)
    rotationM2 = matrix(c(cos(phi.r[2])*cos(theta.r[2]),-sin(phi.r[2]),
                          cos(phi.r[2])*sin(theta.r[2]),sin(phi.r[2])*cos(theta.r[2]),
                          cos(phi.r[2]),sin(phi.r[2])*sin(theta.r[2]),-sin(theta.r[2]),
                          0,cos(theta.r[2])),nrow=3)
    
    fibVec_DiST[rep,1,] = c(-99,-99,-99)
    fibVec_DiST[rep,2,] = c(-99,-99,-99)
    
    phi1_DiST[rep] = -99
    theta1_DiST[rep] = -99
    phi2_DiST[rep] = -99
    theta2_DiST[rep] = -99
    
    fib_true_DiST[rep,1,] = temp$fod1.s%*%rotationM1
    fib_true_DiST[rep,2,] = temp$fod2.s%*%rotationM2
    
    angle1_DiST[rep] = -99
    angle2_DiST[rep] = -99
    sep_DiST[rep] = -99
  }
  
  fit_DiST[[rep]] = out
}

correct_DiST = sum(nfib_DiST==2)/N_rep
over_DiST = sum(nfib_DiST>2)/N_rep
under_DiST = sum(nfib_DiST<2)/N_rep
angle1_DiST_mean = mean(angle1_DiST[nfib_DiST==2])
angle2_DiST_mean = mean(angle2_DiST[nfib_DiST==2])
angle1_DiST_sd = sd(angle1_DiST[nfib_DiST==2])
angle2_DiST_sd = sd(angle2_DiST[nfib_DiST==2])
mean(sep_DiST[nfib_DiST==2])

save(nfib_DiST, angle1_DiST, angle2_DiST, sep_DiST, theta1_DiST, theta2_DiST,
     phi1_DiST, phi2_DiST, fibVec_DiST, fib_true_DiST, fit_DiST,
     correct_DiST, over_DiST, under_DiST,
     file=paste0(path_load,folder_name,c('DiST_summary.RData'),collapse=""))

mat_file = paste0(path_load,folder_name,c('DiST_Matlab_summary.mat'))
writeMat(mat_file, nfib_DiST=nfib_DiST, angle1_DiST=angle1_DiST, 
         angle2_DiST=angle2_DiST, sep_DiST=sep_DiST, theta1_DiST=theta1_DiST, 
         theta2_DiST=theta2_DiST, phi1_DiST=phi1_DiST, phi2_DiST=phi2_DiST,
         fibVec_DiST=fibVec_DiST, fib_true_DiST=fib_true_DiST, 
         correct_DiST=correct_DiST, over_DiST=over_DiST, under_DiST=under_DiST)

# rm(list=ls())