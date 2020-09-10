#ifndef MAIN_H
#define MAIN_H

#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <time.h>

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

/*Creates a folder in current path*/
void create_directory(std::string filename)
{
    char *dir = const_cast<char *>(filename.c_str());
    //mkdir and _mkdir take different arguments, so we need to make sure we are passing them correctly
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    createFolder(dir);
#else
    createFolder(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

/*Returns path to current directory*/
std::string get_current_dir()
{
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);

    for (int i = 0; i < int(current_working_dir.length()); i++)
    {
        if (current_working_dir[i] == '\\')
        {
            current_working_dir.replace(i, 1, slash);
        }
    }
    return current_working_dir;
}

/*Output data to file with the .csv type*/
void writeToFile(std::string filename, int n, double *v, double *u, double *rel_err = 0)
{
    std::ofstream ofile;
    std::string file = filename;

    //Create a new folder in current directory and get path
    create_directory(filename);
    std::string dir = get_current_dir();

    //Append filename with value of n, and give type; i.e filename_n.csv
    file.append("_" + std::to_string(n) + ".csv");

    ofile.open(dir + slash + filename + slash + file, std::fstream::out);

    if (ofile.is_open())
    {
        ofile << "Numeric,";

        if (rel_err != 0)
        {
            ofile << "Analytic,";
            ofile << "Relative Error:" << std::endl;
            for (int i = 0; i < n; i++)
            {
                ofile << v[i] << ",";
                ofile << u[i] << ",";
                ofile << rel_err[i] << std::endl;
            }
        }
        else
        {
            ofile << "Analytic" << std::endl;
            for (int i = 0; i < n; i++)
            {
                ofile << v[i] << ",";
                ofile << u[i] << std::endl;
            }
        }
    }

    ofile.close();
}

/*Writes execution times to file*/
void writeExecTimeToFile(std::string filename, int n, double exec_time)
{
    std::ofstream ofile;
    std::string file = filename;

    exec_time *= 1000; //converting to ms

    create_directory(filename);
    std::string dir = get_current_dir();

    file.append("_" + std::to_string(n) + "_timing" + ".csv");

    std::string fileDir = dir + slash + filename + slash + file;
    //access checks whether or not the file already exists
    if (access(fileDir.c_str(), F_OK) == 0)
    {
        ofile.open(fileDir, std::fstream::app);

        if (ofile.is_open())
        {
            ofile << n << ",";
            ofile << std::fixed << std::setprecision(4) << exec_time << std::endl;
        }
    }
    else
    {
        ofile.open(fileDir, std::fstream::out);

        if (ofile.is_open())
        {
            ofile << "n,";
            ofile << std::fixed << std::setprecision(4) << "ms" << std::endl;
            ofile << n << ",";
            ofile << exec_time << std::endl;
        }
    }
    ofile.close();
}

/*
    This function solves the general equation Av = h^2 *f, using a general algortihm

    a_v --- vector of constants a_1...a_n-1
    b_v --- vector of constants b_1...b_n
    c_v --- vector of constants c_1...c_n-1
    b_tilde = h^2*f(x_i)
    v - unknown
    f - known function
    A - tridiagonal matrix

    n - number of integration points
    h - step size
    */
double *generalSolver(int n, double h, int a, int b, int c)
{
    /*Note: Declaring a_v, b_v and c_v as arrays and filling them with a singular value is technically 
   redundant since the program can only take in one value for each, but since the performance hit is minimal, we will leave it as is.*/
    double *a_v = new double[n - 1];
    double *b_v = new double[n];
    double *c_v = new double[n - 1];

    double *b_tilde = new double[n];
    double *v = new double[n];

    b_v[0] = b;
    b_tilde[0] = h * h * 100 * std::exp(-10 * h);

    for (int i = 0; i < n - 1; i++)
    {
        a_v[i] = a;
        b_v[i + 1] = b;
        c_v[i] = c;
        b_tilde[i + 1] = h * h * 100 * exp(-10 * (i + 2) * h); //i + 2 to compensate for b_tilde[i + 1]
    }

    //Forwards substitution
    for (int i = 1; i < n; i++)
    {
        b_v[i] = b_v[i] - c_v[i - 1] * a_v[i - 1] / b_v[i - 1];
        b_tilde[i] = b_tilde[i] - b_tilde[i - 1] * a_v[i - 1] / b_v[i - 1];
    }

    //Apply backward substitution
    v[n - 1] = b_tilde[n - 1] / b_v[n - 1];
    for (int i = n - 2; i > -1; i--)
    {
        v[i] = (b_tilde[i] - c_v[i] * v[i + 1]) / b_v[i];
    }

    delete[] a_v, b_v, c_v, b_tilde;

    return v;
}

/*
This function solves the general equation Av = h^2 *f for the specific tridiagonal case of -1, 2, -1
For details, see generalSolver()
*/
double *specSolver(int n, double h)
{
    double *b_v = new double[n];
    double *b_tilde = new double[n];
    double *v = new double[n];

    b_v[0] = 2;
    b_tilde[0] = h * h * 100 * std::exp(-10 * h);

    for (int i = 1; i < n; i++)
    {
        b_v[i] = (i + 2.0) / (i + 1.0);                                                  //analytical expression of b_v
        b_tilde[i] = h * h * 100 * exp(-10 * (i + 1) * h) + b_tilde[i - 1] / b_v[i - 1]; //immediately calculating b_tilde with forward substitution
    }

    v[n - 1] = b_tilde[n - 1] / b_v[n - 1]; //Solving the final row

    //Apply backward substitution
    for (int i = n - 2; i > -1; i--)
    {
        v[i] = (b_tilde[i] + v[i + 1]) / b_v[i];
    }

    delete[] b_v, b_tilde;

    return v;
}

/*
Using Armadillo for LU decomp. We have A*v = b_tilde
Since A = L*U, we define U*v = w, which means L*w = b_tilde
Then we solve first w and then v.
*/
double *lusolver(int n, double h, int a, int b, int c)
{
    //using armadillo for the matrix handling
    arma::vec b_tilde(n);
    arma::vec w(n);
    arma::vec v(n);
    arma::mat A(n, n);
    double *v_pointer = new double[n]; //for compatibility with writeToFile

    //initializing the matrix A with -1 2 -1
    A.fill(0.0);
    A(0, 0) = b;
    /*There is a better way of doing this, using armadillo A.diag()*/
    for (int i = 1; i < n; i++)
    {
        A(i, i) = b;
        A(i - 1, i) = c;
        A(i, i - 1) = a;
    }

    //LU decomp
    arma::mat L, U;
    arma::lu(L, U, A);

    //solving w
    for (int i = 0; i < n; i++)
    {
        b_tilde(i) = h * h * 100 * std::exp(-10 * h * (i + 1));
        w(i) = b_tilde(i);
        for (int j = 0; j < i; j++)
        {
            w(i) -= L(i, j) * w(j);
        }
    }

    //solving v
    for (int i = n - 1; i > -1; i--)
    {
        v(i) = w(i);
        for (int j = n - 1; j > i; j--)
        {
            v(i) -= U(i, j) * v(j);
        }
        v(i) /= U(i, i);
        v_pointer[i] = v(i);
    }

    return v_pointer;
}

/*Returns the analytical solution of equation u, with n integration points*/
double *analyticalSolution(int n, double h)
{
    double *u = new double[n];
    for (int i = 0; i < n; i++)
    {
        u[i] = 1 - (1 - std::exp(-10)) * (i + 1) * h - std::exp(-10 * (i + 1) * h);
    }
    return u;
}

/*Times and outputs the execution of the algorithms*/
void time_and_write(double *(*solver)(int, double, int, int, int), int n, int a, int b, int c, std::string task)
{
    std::clock_t start, finish;
    double h = 1.0 / (n + 1);
    double *v;

    if (a == 0 || b == 0 || c == 0)
    {
        start = std::clock();
        v = specSolver(n, h);
        finish = std::clock();
    }
    else
    {
        start = std::clock();
        v = solver(n, h, a, b, c);
        finish = std::clock();
    }

    double *u = analyticalSolution(n, h);

    double execution_time = double(finish - start) / double(CLOCKS_PER_SEC);
    std::cout << "Execution time for n = " << n << " is " << std::fixed << std::setprecision(4)
              << execution_time * 1000 << "ms" << std::endl;

    writeToFile(task, n, v, u);
    writeExecTimeToFile(task, n, execution_time);

    delete[] v, u;
}

#endif