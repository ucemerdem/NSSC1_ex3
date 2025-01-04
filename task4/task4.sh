#!/bin/bash
# Number of cores
#SBATCH -c 1
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS 
export OMP_NUM_THREADS=1
# You can start several programs with one script file/submission
./solver_static static 2048 100000
./solver_static1 static1 2048 100000
./solver_dynamic  dynamic 2048 100000



