#!/usr/bin/env python3
"""Assignment 1: Parallelisation and Improvement of code to determine Pi."""
import time
from mpi4py import MPI
from mpmath import quad
start_time = time.time()
comm = MPI.COMM_WORLD

nproc = comm.Get_size()
# The first processor is leader, so one fewer available to be a worker
nworkers = nproc - 1

# samples
N = 16
DELTA = 1.0 / N

# integral
I = 0.0

def integrand(x):
    """Define Integrand"""
    return 4.0 / (1.0 + x*x)

if comm.Get_rank() == 0:

  # Leader: choose points to sample function, send to workers and
  # collect their contributions. Also calculate a sub-set of points.

    for i in range(0,N):

    # decide which rank evaluates this point
        j = i % nproc

    # mid-point rule
        z = (i+0.5) * DELTA

        if j == 0:
      # so do this locally using the leader machine
            y = quad(integrand, [i * DELTA, (i+1) * DELTA])
        else:
      # communicate to a worker
            comm.send(z, dest=j)
            y = comm.recv(source=j)

        I += y

  # Shut down the workers
    for i in range(1, nproc):
        comm.send(-1.0, dest=i)

    print(f"Integral % {I}")

else:

  # Worker: waiting for something to happen, then stop if sent message
  # outside the integral limits

    while True:

        z = comm.recv(source=0)

        if z < 0.0:
      # stop the worker
            break
        comm.send(quad(integrand, [i * DELTA, (i+1) * DELTA]), dest=0)

end_time = time.time()
print("The runtime was", f"{end_time - start_time}","seconds")
