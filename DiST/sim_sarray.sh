#!/bin/bash -l
# NOTE the -l flag!

###############################################################################
##
## NOTES:
##
## Submit as:
##
##    sarray ./sim_sarray.sh
##
##    or, for array jobs:
##
##    sbatch ./sim_sarray.sh
##
## (1) When specifying --range as a range it must start from a positive
##     integer e.g.,
##       SARRAY --range=0-9
##     is not allowed.
##
## (2) Negative numbers are not allowed in --range
##     e.g.,
##      SARRAY --range=-5,-4,-3,-2,-1,0,1,2,3,4,5
##     is not allowed.
##
## (3) Zero can be included if specified separately.
##    e.g.,
##       SARRAY --range=0,1-9
##     is allowed.
##
## (4) Ranges can be combined with specified job numbers.
##    e.g.,
##       SARRAY --range=0,1-4,6-10,50-100
##     is allowed.
##
###############################################################################

# Load R module:
module load R
export R_LIBS=/home/rkwwong/R/x86_64-pc-linux-gnu-library/3.0:$R_LIBS

# Name of the job - You'll probably want to customize this.
#SBATCH --job-name=tt
#SBATCH --mem-per-cpu=1500
#SBATCH --exclude=c0-[10-11]
#SBATCH --cpus-per-task=4

# Array job specifications:
#SARRAY --range=1-200

# Email notifications (optional), type=BEGIN, END, FAIL, ALL
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=rkwwong@ucdavis.edu

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o log/test_%j.out
#SBATCH -e log/test_%j.err

# Execute each of the jobs with a different index (the R script will then process
# this to do something different for each index):
srun --mem-per-cpu=1500 --exclude=c0-[10-11] --cpus-per-task=4 R --vanilla --no-save --args ${SLURM_ARRAYID} < test.R





