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
#SBATCH --ntasks=8
# Runtime (hard, HH:MM:SS)
#SBATCH --time=00:30:00
# Job name
#SBATCH --job-name=assignment_1_naive
# Output file
#SBATCH --output=assignment_1_naive-%j.out
# choose which version to load

module purge
module load openmpi/gcc-8.5.0/4.1.1
/opt/software/scripts/job_prologue.sh


pylint --extension-pkg-allow-list=mpi4py.MPI assignment_1_naive
# Modify the line below to run your program
mpirun -np $SLURM_NPROCS ./assignment_1_naive

/opt/software/scripts/job_epilogue.sh
