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
./mcint SINX -1.570796 1.570796 1e7
./mcint COS2XINV -1.570796 1.570796 1e7
./mcint X4M5 -1.570796 1.570796 1e7
