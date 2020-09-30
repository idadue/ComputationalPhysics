# Project 2: Eigenvalue problems

## The report

Link to the report when it is ready

## How to run the code
The code can be compiled by e.g.: 

```
g++ -o main.out main.cpp jacobisolver.cpp -larmadillo -llapack -lblas
```

The program reads data from the .txt file inData.txt.

The format of the file is as such:
```
n rho potential omega folder_name arma
```

The program supports a number of lines after each other, and ignores all text before it reaches the special character `\`. In place of `arma`, one can use any charachter instead. In that case, the program will not compute the eigenvalues and eigenvectors using armadillo functions. If `folder_name` is given as simply `.`, then the program will not write the results to file and for any other string, results will be written in /results/folder_name.

### 1. Buckling Beam problem:
  To solve for the buckling beam problem, use potential = 0 in inData.txt
  
  Example:
```
    \100 2 0 0 bbeam arma
```
   
  This will run the program once for n = 100, rho = 2, saving results to file and also compute with `armadillo`.

### 2. Quantum dots with one electron:
  To solve for the quantum dots with one electron problem, use potential = 1 in inData.txt
  
  Example:
```
    \100 2 1 0 . arma
```
   
  This will run the program once for n = 100, rho = 2, not save results to file and also compute with `armadillo`.


### 3. Quantum dots with two electrons: 
  To solve for the quantum dots with one electron problem, use potential = 1 in inData.txt
  
  Example:
```
    \100 10 1 0.5 . .
```
   
  This will run the program once for n = 100, rho = 10, not save results to file and not compute with armadillo.
 
An example of stacking is show here:
inData.txt
```
Anything written her will not be run..
\
100 2 0 0 bbeam arma
100 2 0 0 bbeam arma
100 2 0 0 bbeam arma
100 2 0 0 bbeam arma
100 2 0 0 bbeam arma

100 5 1 0 . arma
100 6 1 0 . arma
100 7 1 0 . arma
100 8 1 0 . arma

100 20 1 0.01 . .
100 20 1 0.1 . .
100 20 1 0.5 . .
100 20 1 1 . .

end inData.txt
```
The program will run through all these scenarios.
