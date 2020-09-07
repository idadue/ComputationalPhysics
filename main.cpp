#include <iostream>
#include <fstream>
#include <string>
#include "time.h"
#include <iomanip>
#include <cmath>
#include <algorithm>

/*
* This part checks during pre proccessing what libraries are included in the users system. If the user is using windows, they should have
* direct.h located on their system, and if the user has a unix based system, they should have unistd.h and sys/stat.h 
* on their system. We only want to include libraries that the user already has, and ignore the others, such that we do not get any compile errors.
*/
#if __has_include(<direct.h>)
#include <direct.h>
#define GetCurrentDir _getcwd
#define createFolder _mkdir

#else
#include <unistd.h>
#include <sys/stat.h>
#define createFolder mkdir
#define GetCurrentDir getcwd
#endif

//Here we define what value the string "slash" should have, dependent on what system the user is using
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
std::string slash = "//";
#else
std::string slash = "/";
#endif

void create_directory(std::string filename) {
    char* dir = const_cast<char*>(filename.c_str());    //convert string to char*
    //mkdir and _mkdir take different arguments, so we need to make sure we are passing them correctly
    #if defined(__CYGWIN__)
    createFolder(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #else
    createFolder(dir);
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

double* generalSolver(int n, double h, int a, int b, int c) {
    /*
    This function solves the general equation Av = h^2 *f

    a_v --- vector of constants a_1...a_n-1
    b_v --- vector of constants b_1...b_n
    c_v --- vector of constants c_1...c_n-1
    b_tilde = h²*f(x_i)

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

void writeToFile(std::string filename, int n, double* v, double* u, double* rel_err = 0) {
    //Tries to create a folder in current directory and create and write to files in that new folder

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

//TODO Find also the precise number of floating point operations needed to solve the above equations.

//Begin main program
int main(int argc, char* argv[]) {
    // A*v = h^2 *f_i = b_tilde
    //A needs a, b and c values. 
    //if constant, they could be hard coded
    // v is what we will calculate
    // h is found from n, which are our integration points
    // b_tilde is found from f, which we will calculate

    //Assumed as integers in this case, but we will treat like arrays in the implementation
    //If we wanted to use n different integers we would have to rewrite some of the code, and it would be better to
    //read such an array from a file
    //We could also just have the program check if it is given an array of numbers or just a single
    //int, and then do different actions based on the situation
    int a;
    int b;
    int c;

    if (argc < 4) {
        //If no commands are specified, exit program
        //Could instead provide a standard set of values
        std::cout << "Bad Usage: " << argv[0] << " is missing arguments" << std::endl;
        exit(1);
    }
    else {
        a = std::atoi(argv[1]);
        b = std::atoi(argv[2]);
        c = std::atoi(argv[3]);
        char task = ' ';

        std::clock_t start, finish;

        std::cout << "Valid tasks are b, c, d, or e, f for custom and 0 to exit. " << std::endl;
        while (task != '0') {
            std::cout << "Insert task: ";
            std::cin >> task;

            switch (task) {
            //TODO declare a global n array for all tasks?
            case 'b': {
                int n[3] = { 10, 100, 1000 };
                for (int i = 0; i < 3; i++) {
                    double h = 1.0 / (n[i] + 1);

                    //Timing the execution of the algorithm
                    start = std::clock();
                    double* v = generalSolver(n[i], h, a, b, c);
                    finish = std::clock();

                    double* u = analyticalSolution(n[i], h);

                    double execution_time = double(finish - start) / double(CLOCKS_PER_SEC);
                    std::cout << "Execution time for n = " << n[i] << " is " << std::fixed << std::setprecision(4) 
                        << execution_time*1000 << "ms" << std::endl;

                    writeToFile("task_b", n[i], v, u);

                    delete[] v, u;
                    v, u = NULL;
                }
                std::cout << "Task b has been completed \n" << std::endl;
                break;
            }
            case 'c': {
                //TODO: Create a specialized algo
                //Similar structure as task b
                int n[5] = { 10, 100, 1000, 10000, 1000000 };
                std::cout << "Not yet implemented" << std::endl;
                exit(1);
            }
            case 'd': {
                //TODO not finished yet
                int n = int(pow(10, 7));
                double h = 1.0 / (n + 1);

                //Can use specialized algorithm instead
                double* v = generalSolver(n, h, a, b, c);
                double* u = analyticalSolution(n, h);

                //The realtive error
                double* rel_err = new double[n];
                for (int i = 0; i < n; i++) {
                    rel_err[i] = std::log10(std::abs((v[i] - u[i]) / u[i]));
                }
                std::cout << "Writing to file, this might take a while." << std::endl;
                writeToFile("task_c", n, v, u, rel_err);

                std::cout << "Not completely implemented yet" << std::endl;

                delete[] v, u;
                v, u = NULL;
                break;
            }
            case 'e': {
                //Use armadillo or make our own implementation of LU decomp. ?
                int n[3] = { 10, 100, 1000 };
                std::cout << "Not yet implemented" << std::endl;
                exit(1);
            }
            case 'f': {
                int n; std::string filename;
                std::cout << "Enter value for n, filename: ";
                std::cin >> n >> filename;
                std::cout << "\n";
                std::cout << "You entered: " << n << ", " << filename << std::endl;

                double h = 1.0 / (n + 1);

                //Timing the execution of the algorithm
                start = std::clock();
                double* v = generalSolver(n, h, a, b, c);
                finish = std::clock();

                double* u = analyticalSolution(n, h);

                double execution_time = double(finish - start) / double(CLOCKS_PER_SEC);
                std::cout << "Execution time for n = " << n << " is " << std::fixed << std::setprecision(4)
                    << execution_time * 1000 << "ms" << std::endl;

                std::cout << "Writing data to file..." << std::endl;
                writeToFile(filename, n, v, u);

                delete[] v, u;
                v, u = NULL;

                std::cout << "Task completed! \n" << std::endl;
                break;
            }
            default: {
                if (task == '0') {
                    std::cout << "Exiting..." << std::endl;
                }
                else {
                    std::cout << "Error, must insert valid character!" << std::endl;
                    break;
                }
            }

            }

        }
    }



    return 0;
}