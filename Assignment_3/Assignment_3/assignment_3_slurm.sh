#!/bin/bash
# Propagate environment variables to the compute node
#SBATCH --export=ALL
# Run in the standard partition (queue)
#SBATCH --partition=teaching
# Specify project account
#SBATCH --account=teaching
# Distribute processes in round-robin fashion
#SBATCH --distribution=block:block
# No of cores required (max. of 16, 4GB RAM per core)
#SBATCH --ntasks=16
# Runtime (hard, HH:MM:SS)
#SBATCH --time=00:30:00
# Job name
#SBATCH --job-name=assignment_3_code.py
# Output file
#SBATCH --output=assignment_3_code.out
# choose which version to load

module purge
module load openmpi/gcc-8.5.0/4.1.1
module load miniconda/3.12.8
/opt/software/scripts/job_prologue.sh

pylint --extension-pkg-allow-list=mpi4py.MPI assignment_3_code.py
# Modify the line below to run your program
mpirun -np 1 ./assignment_3_code.py

#mpirun -np 2 ./assignment_3_code.py

#mpirun -np 4 ./assignment_3_code.py

#mpirun -np 8 ./assignment_3_code.py

#mpirun -np 16 ./assignment_3_code.py
# Output
/opt/software/scripts/job_epilogue.sh
