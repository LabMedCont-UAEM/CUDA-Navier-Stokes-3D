# CUDA-Navier-Stokes-3D
A complete Navier-Stokes's equation solver using diferrent geometry mesh, parallelization is available by CUDA.

You must to have the PGI's group compiler to run CUDA Fortran codes.

Compiling:
$ pgfortran CUDA_NS3D_Cart.f90 -Mcuda=fastmath -o CUDANS3D
