# Project 4: The Ising model

To run the code using MPI you need to have MPI installed on your computer. 
On ubuntu you should be able to run ```sudo apt install libopenmpi-dev```
  
Once you have MPI, run 
  ```mpiCC --showme```
in the terminal to see where on your computer the include and lib folders are located. Replace the paths from your system with the ones already in the makefile in /src/MPI.

Then just run ```make``` with whatever values are desired.

Running the code in /no_mpi is as simple as just running ```make``` in terminal.
