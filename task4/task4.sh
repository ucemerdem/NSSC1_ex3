#!/bin/bash

threads=(1 2 4 8 10 20 40)
cores=(1 1 2 4 5 10 20)

for i in "${!threads[@]}"; do
    job_script="job_${threads[$i]}threads_${cores[$i]}cores.sh"

    cat << EOF > $job_script
#!/bin/bash
#SBATCH -c ${cores[$i]}
#SBATCH --time=01:00:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS
export OMP_NUM_THREADS=${threads[$i]}
# Run your programs
./solver_static static 2048 100000
./solver_static1 static1 2048 100000
./solver_dynamic dynamic 2048 100000
EOF

    sbatch $job_script

done