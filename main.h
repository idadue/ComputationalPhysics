#ifndef MAIN_H
#define MAIN_H

#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include <armadillo>

/*
* This part checks during pre proccessing what system the user is using. If the user is using windows, they should have
* direct.h and io.h located on their system, and if the user has a unix based system, they should have unistd.h and sys/stat.h
* on their system. We only want to include libraries that the user already has, and ignore the others, such that we do not get any compile errors.
*/

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#include <io.h>
#define access _access_s
#define GetCurrentDir _getcwd
#define createFolder _mkdir
#define F_OK 0
std::string slash = "//";
#else
#include <unistd.h>
#include <sys/stat.h>
#define createFolder mkdir
#define GetCurrentDir getcwd
std::string slash = "/";
#endif

//Directory finders/editors
void create_directory(std::string filename) {
    char* dir = const_cast<char*>(filename.c_str());    //convert string to char*
    //mkdir and _mkdir take different arguments, so we need to make sure we are passing them correctly
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        createFolder(dir);
    #else
        createFolder(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif
}

std::string get_current_dir() {
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);

    for (int i = 0; i < int(current_working_dir.length()); i++) {
        if (current_working_dir[i] == '\\') {
            current_working_dir.replace(i, 1, slash);
        }
    }
    return current_working_dir;
}

//Output functions
void writeToFile(std::string filename, int n, double* v, double* u, double* rel_err = 0) {
    //Prepare output to file
    std::ofstream ofile;

    //Create new file
    std::string file = filename;

    //Create a new folder in current directory
    create_directory(filename);
    std::string dir = get_current_dir();

    //Append filename with value of n, and give type; i.e filename_n.csv
    file.append("_" + std::to_string(n) + ".csv");

    ofile.open(dir + slash + filename + slash + file, std::fstream::out);
    //Write to file here
    if (ofile.is_open()) {
        ofile << "Numeric,";

        if (rel_err != 0) {
            ofile << "Analytic,";
            ofile << "Relative Error:" << std::endl;
            for (int i = 0; i < n; i++) {
                ofile << v[i] << ",";
                ofile << u[i] << ",";
                ofile << rel_err[i] << std::endl;
            }
        }
        else {
            ofile << "Analytic" << std::endl;
            for (int i = 0; i < n; i++) {
                ofile << v[i] << ",";
                ofile << u[i] << std::endl;
            }
        }
    }

    ofile.close();
}

void writeExecTimeToFile(std::string filename, int n, double exec_time) {
    //Writes execution times to file
    std::ofstream ofile;
    std::string file = filename;

    exec_time *= 1000; //converting to ms

    create_directory(filename);
    std::string dir = get_current_dir();

    file.append("_" + std::to_string(n) + "_timing" + ".csv");

    std::string fileDir = dir + slash + filename + slash + file;
    if (access(fileDir.c_str(),F_OK) == 0) {
        ofile.open(fileDir, std::fstream::app);

        if (ofile.is_open()) {
            ofile << n << ",";
            ofile << std::fixed << std::setprecision(4) << exec_time << std::endl;
        }
    }
    else {
        ofile.open(fileDir, std::fstream::out);

        if (ofile.is_open()) {
            ofile << "n,";
            ofile << std::fixed << std::setprecision(4) << "ms" << std::endl;
            ofile << n << ",";
            ofile << exec_time << std::endl;
        }
    }
    ofile.close();
}

//Solver functions
double* generalSolver(int n, double h, int a, int b, int c) {
    /*
    This function solves the general equation Av = h^2 *f

    a_v --- vector of constants a_1...a_n-1
    b_v --- vector of constants b_1...b_n
    c_v --- vector of constants c_1...c_n-1
    b_tilde = h�*f(x_i)

    n - number of integration points
    h - step size
    */


    double* a_v = new double[n - 1];
    double* b_v = new double[n];
    double* c_v = new double[n - 1];

    double* b_tilde = new double[n];
    double* v = new double[n]; //the vector we're trying to solve, final answer

    b_v[0] = b;
    b_tilde[0] = h * h * 100 * std::exp(-10 * h);

    for (int i = 0; i < n - 1; i++) {
        a_v[i] = a;
        b_v[i + 1] = b;
        c_v[i] = c;
        b_tilde[i + 1] = h * h * 100 * exp(-10 * (i + 2) * h);  //i + 2 to compensate for b_tilde[i + 1]
    }

    for (int i = 1; i < n; i++) {
        b_v[i] = b_v[i] - c_v[i - 1] * a_v[i - 1] / b_v[i - 1];
        b_tilde[i] = b_tilde[i] - b_tilde[i - 1] * a_v[i - 1] / b_v[i - 1];
    }

    v[n - 1] = b_tilde[n - 1] / b_v[n - 1]; //Solving the final row

    //Apply backward substitution
    for (int i = n - 2; i > -1; i--) {
        v[i] = (b_tilde[i] - c_v[i] * v[i + 1]) / b_v[i];
    } //Solving all other rows

    delete[] a_v, b_v, c_v, b_tilde;
    a_v, b_v, c_v, b_tilde = NULL;

    return v;
}

double* specSolver(int n, double h) {
    /*
    This function solves the general equation Av = h^2 *f for the specific tridiagonal case of -1, 2, -1

    b_v --- vector of constants b_1...b_n, given by b_i = (i+1)/i
    b_tilde = h²*f(x_i)

    n - number of integration points
    h - step size
    */


    double* b_v = new double[n];
    double* b_tilde = new double[n];
    double* v = new double[n]; //the vector we're trying to solve, final answer

    b_v[0] = 2;
    b_tilde[0] = h * h * 100 * std::exp(-10 * h);

    for (int i = 1; i < n; i++) {
      b_v[i] = (i + 2.0)/(i + 1.0); //analytical expression of b_v
      b_tilde[i] = h * h * 100 * exp(-10 * (i + 1) * h) + b_tilde[i - 1] / b_v[i - 1]; //immediately calculating b_tilde with forward substitution
    }

    v[n - 1] = b_tilde[n - 1] / b_v[n - 1]; //Solving the final row

    //Apply backward substitution
    for (int i = n - 2; i > -1; i--) {
        v[i] = (b_tilde[i] + v[i + 1]) / b_v[i];
    } //Solving all other rows

    delete[] b_v, b_tilde;
    b_v, b_tilde = NULL;

    return v;
}

double* lusolver(int n, double h, int a, int b, int c) {

    /*
    Using Armadillo for LU decomp. We have A*v = b_tilde
    Since A = L*U, we define U*v = w, which means L*w = b_tilde
    Then we solve first w and then v.
    */


    //using armadillo for the matrix handling
    arma::vec b_tilde(n);
    arma::vec w(n);
    arma::vec v(n);
    arma::mat A(n,n);
    double* v_pointer = new double[n]; //for compatibility with writeToFile

    //initializing the matrix A with -1 2 -1
    A.fill(0.0);
    A(0,0) = b;
    for ( int i = 1; i < n; i++){
      A(i,i) = b;
      A(i-1,i) = c;
      A(i,i-1) = a;
    }

    //LU decomp
    arma::mat L, U;
    arma::lu(L, U, A);

    /*
    printing for testing purposes
    L.print("L= ");
    U.print("U= ");
    A.print("A =");
    (A-L*U).print("LU test");
    */

    //solving w
    for ( int i = 0; i < n; i++){
      b_tilde(i) = h * h * 100 * std::exp(-10 * h * (i+1));
      w(i) = b_tilde(i);
      for ( int j = 0; j < i ; j++){
        w(i) -= L(i,j)*w(j);
      }
    }

    //solving v
    for ( int i = n - 1; i > -1; i--){
      v(i) = w(i);
      for ( int j = n - 1; j > i ; j--){
        v(i) -= U(i,j)*v(j);
      }
      v(i) /= U(i,i);
      v_pointer[i] = v(i);
    }

    return v_pointer;
}

//Analytical solution
double* analyticalSolution(int n, double h) {
    /*
    Returns the analytical solution of equation u, with n integration points
    */

    double* u = new double[n];
    for (int i = 0; i < n; i++) {
        u[i] = 1 - (1 - std::exp(-10)) * (i + 1) * h - std::exp(-10 * (i + 1) * h);
    }
    return u;
}


#endif
