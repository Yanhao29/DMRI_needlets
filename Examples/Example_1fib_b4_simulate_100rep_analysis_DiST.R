library(R.matlab)

path_load = '/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/'
folder_name = '1fib_lmax8_b4_ratio10_n41_sig0.05/'
file_name = '1fib_lmax8_b4_ratio10_n41_sig0.05_rep'

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

theta1_DiST = rep(0,N_rep)
phi1_DiST = rep(0,N_rep)

fibVec_DiST = matrix(0,nrow = N_rep, ncol = 3)
fib_true_DiST = matrix(0,nrow = N_rep, ncol = 3)

fit_DiST = rep(list(),N_rep)

for(rep in 1:N_rep){
  ## load matlab simulation
  temp = readMat(paste0(path_load,folder_name,file_name,toString(rep),'.mat'))
  grad.mat_temp <- t(temp$pos.sampling)
  design_temp <- generate.design.noS0(grad.mat_temp)
  
  print(rep)
  
  out <- try(biased.mle.cs2.iso(dwis=c(temp$DWI), design=design_temp,
                                sigma=c(temp$sig), S0=1, b=c(temp$b), maxJ=opt$maxJ,
                                basis.depth=opt$basis.depth,
                                basis.random=opt$basis.random,
                                def.alphas=opt$def.alphas, betas0=NULL,
                                Ninner=opt$Ninner, Nmiddle.max=opt$Nmiddle.max,
                                display=F, pre.thres=opt$pre.thres,
                                alphas.const=opt$alphas.const))
  nfib_DiST[rep] = out$n.fiber
  
  if(out$n.fiber == 1){
    fib_temp = out$vs
    theta.r = temp$theta.r
    phi.r = temp$phi.r
  
    fibVec_DiST[rep,] = out$vs
    phi1_DiST[rep] = temp$phi.r
    theta1_DiST[rep] = temp$theta.r
    rotationM = matrix(c(cos(phi.r)*cos(theta.r),-sin(phi.r),cos(phi.r)*sin(theta.r),sin(phi.r)*cos(theta.r),cos(phi.r),sin(phi.r)*sin(theta.r),-sin(theta.r),0,cos(theta.r)),nrow=3)
    
    fib_true = temp$fod1.s%*%rotationM
    fib_true_DiST[rep,] = fib_true
    
    angle1_DiST[rep] = min(acos(sum(fib_temp*fib_true)/(sqrt(sum(fib_temp*fib_temp))*sqrt(sum(fib_true*fib_true)))),
                           abs(pi-acos(sum(fib_temp*fib_true)/(sqrt(sum(fib_temp*fib_temp))*sqrt(sum(fib_true*fib_true))))))*180/pi
  } else {
    theta.r = temp$theta.r
    phi.r = temp$phi.r
    rotationM = matrix(c(cos(phi.r)*cos(theta.r),-sin(phi.r),cos(phi.r)*sin(theta.r),sin(phi.r)*cos(theta.r),cos(phi.r),sin(phi.r)*sin(theta.r),-sin(theta.r),0,cos(theta.r)),nrow=3)
    
    fibVec_DiST[rep,] = c(-99,-99,-99)
    phi1_DiST[rep] = -99
    theta1_DiST[rep] = -99
    fib_true_DiST[rep,] = temp$fod1.s%*%rotationM
    angle1_DiST[rep] = -99
  }
  
fit_DiST[[rep]] = out
}

correct_DiST = sum(nfib_DiST==1)/N_rep
over_DiST = sum(nfib_DiST>1)/N_rep
under_DiST = sum(nfib_DiST<1)/N_rep
angle1_DiST_mean = mean(angle1_DiST[nfib_DiST==1])

save(nfib_DiST, angle1_DiST, theta1_DiST, phi1_DiST, fibVec_DiST, fib_true_DiST, fit_DiST,
     correct_DiST, over_DiST, under_DiST,
     file=paste0(path_load,folder_name,c('DiST_summary.RData'),collapse=""))

mat_file = paste0(path_load,folder_name,c('DiST_Matlab_summary.mat'))
writeMat(mat_file, nfib_DiST=nfib_DiST, angle1_DiST=angle1_DiST, theta1_DiST=theta1_DiST,
         phi1_DiST=phi1_DiST, fibVec_DiST=fibVec_DiST, fib_true_DiST=fib_true_DiST,
         correct_DiST=correct_DiST, over_DiST=over_DiST, under_DiST=under_DiST)

