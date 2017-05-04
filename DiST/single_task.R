################################################################################
#### Simulation for testing smoothing methods over new prelimilary fittings ####
################################################################################
# for method 30
# need ncpu

################
#### Source ####
################
source("dwi_fit.R")

#######################
#### Simulate data ####
#######################
source("sim-curve-new.R")

####################################
#### fit voxel level estimation ####
####################################
pre <- v.est(dwi.obs=dwi.obs, sigma=sigma, S0=S0const, b=1, grad.mat=grad.mat,
             braingrid=braingrid, cpus=ncpu, opt=list())

save(pre, dwi.obs, tensor.true, n.layer, braindim, braingrid, n.fiber, v0.true,
     file=paste(c("/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/fit.RData"),collapse=""))

