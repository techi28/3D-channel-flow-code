# 3D-channel-flow-code
This code is developed by Miguel and Xin (LHSV), including a test case of 3D channel flow.
It is written in FORTRAN90. It solves the Navier-Stokes equations on a collocated grid and combines the finite volume method with the sigma-transformation method. It is written in sequential or in parallel with the MPI library that can run on multiple processors such as high performance computer.

'src' is the source file.
'input' is the input file for prepocess like generating mesh and divide domain for HPC. 
'output' is the result file.

*********
