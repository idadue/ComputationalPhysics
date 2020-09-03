#include <iostream>
#include <fstream>
#include <string>
#include "time.h"

//include setw, setprecision, etc
#include <iomanip>
#include <cmath> //needed for exp


//begin main program
int main(int argc, char *argv[]){

    //Input from command line
    // A*u = h^2 b_tilde
    //A needs a, b and c values. Could be vector of constants, but need only be single constant for this project?
    //for 1b we need it to be general, so should be vectors - G
    //if constant, they could be hard coded
    // u is what we will calculate
    // h is found from n, which are our integration points
    // b is found from f, which we will calculate

    int n;
    std::string filename;
    int a;
    int b;
    int c;

    //Read in the name of the output file, and the exponent n
    if (argc < 6) {
        //If no commands are specified, exit program
        //Could instead provide a standard set of values
        std::cout << "Bad Usage: " << argv[0] << " is missing arguments" << std::endl;
        exit(1);
    }
    else {
        filename = argv[1]; // first command line argument after name of program
        n = std::atoi(argv[2]);
        a = std::atoi(argv[3]);
        b = std::atoi(argv[4]);
        c = std::atoi(argv[5]);
    }


    //Prepare output to file
    std::ofstream ofile;
    //Create new file
    std::string file = filename;

    //Append filename with value of n, and give type; i.e filename_n.csv
    file.append("_" + std::to_string(n) + ".csv");

    //Calculation starts here

    //Time the execution of the algorithm
    std::clock_t start, finish;
    start = std::clock();

    //Vector versions of the input values from command line
    double * v_a = new double [n-1];
    double * v_b = new double [n];
    double * v_c = new double [n-1];
    double h = 1.0 / (n + 1);

    for (int i = 0; i < n-1; i++) {
      v_a[i] = a;
      v_c[i] = c;
    };
    //Vector value of the right side of the equation, h²*f(x_i)
    double * v_g = new double[n];
    for ( int i = 0; i < n; i++) {
      v_g[i] = h*h*100*exp(-10*(i+1)*h);
      v_b[i] = b;
    };

    double * B = new double[n]; //Diagonal values after forward substitution
    double * G = new double[n]; //b_tilde values after forward substitution
    double * u = new double[n]; //the vector we're trying to solve, final answer

    G[0] = v_g[0];
    B[0] = v_b[0]; //Forward substitution leaves the first element of b and g unchanged

    for ( int i = 1; i < n; i++) {
      B[i] = v_b[i] - v_c[i-1]*v_a[i-1]/B[i-1];
      G[i] = v_g[i] - G[i-1]*v_a[i-1]/B[i-1];
    } //Performing forward substitution

    u[n-1] = G[n-1]/B[n-1]; //Solving the final row

    for ( int i = n-2; i > -1; i--) {
      u[i] = (G[i] - v_c[i]*u[i+1])/B[i];
    } //Solving all other rows

    finish = std::clock();
    std::cout << "Execution time: " << float(finish - start) / CLOCKS_PER_SEC << " seconds" << std::endl;
    //Calculation ends here

    //Open and close file
    ofile.open(file);
    //Write to file here
    //can format the output data later
    ofile << "Numeric:,";
    ofile << "Analytic:,";
    ofile << "Relative Error:" << std::endl;

    double * v = new double[n];
    double * relative_err = new double[n];

    for ( int i = 0; i < n; i++) {
      v[i] = 1 - (1 - std::exp(-10))*(i + 1)*h - std::exp(-10*(i + 1)*h);
      relative_err[i] = ((v[i] - u[i]) / v[i]);
      ofile << u[i] << ",";
      ofile << v[i] << ",";
      ofile <<  relative_err[i] << std::endl;
    }//printing current results compared with the analytic solution

    ofile.close();


    delete[] v_a, v_b, v_c, v_g;
    delete[] B, G, u;
    delete[] relative_err;
    v_a, v_b, v_c, v_g = NULL;
    B, G, u = NULL;
    relative_err = NULL;

    return 0;
}
