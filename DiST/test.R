
#============================== Setup for running on Gauss... ==============================#


"ppaste" <- function(...){paste(...,sep="")}

args <- commandArgs(TRUE)

cat(ppaste("Command-line arguments:\n"))
print(args)

####
# sim_zero ==> Lowest possible dataset number
# sim_start ==> Lowest dataset number to be analyzed by this particular batch job
###

###################
sim_start <- 0
#sim_zero <- 1
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
  sinkit <- FALSE
} else {
  # SLURM can use either 0- or 1-indexing...
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
  set.seed(762*(sim_num-1) + 121231)
}

#i <- sim_num # use i rather than sim_num
sinkfile <- paste("log/output_progress_",sim_num,".txt",sep="")

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

#============================== Run the simulation study ==============================#

if (sinkit){
  cat(paste("Sinking output to: ",sinkfile,"\n",sep=""))
  sink(sinkfile)
}

# Load dataset sim_num...
# load(paste("Data_",sim_num,".RData",sep=""))

## Do lots of cool stuff analyzing dataset sim_num..
## Do not use sim_num

ncpu <- 4 # need to match with sim_sarray.sh
source("single_task.R")


# Save dataset sim_num...
save(pre, dwi.obs, tensor.true, n.layer, braindim, braingrid, n.fiber, v0.true,
     file=paste(c("sim_results/fit-",sim_num,".RData"),collapse=""))


