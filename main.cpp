#include <iostream>
#include <fstream>
#include <string>
#include "time.h"
#include <iomanip>
#include <cmath> 
#include <direct.h> // mkdir
#include <filesystem>

//TODO Find also the precise number of floating point operations needed to solve the above equations.

double* generalSolver(int n, double h, int a, int b, int c) {
    //Function that solves the general equation on
    //the form Av = b

    //a_v --- vector of constants a_1...a_n-1
    //b_v --- vector of constants b_1...b_n
    //c_v --- vector of constants c_1...c_n-1
    //b_tilde = h²*f(x_i)


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

    //Get current path and convert to string
    std::filesystem::path path = std::filesystem::current_path(); //C++17
    std::string current_path(path.u8string());

    //Append filename with value of n, and give type; i.e filename_n.csv
    file.append("_" + std::to_string(n) + ".csv");

    //Creates new folder in directory
    //might not work on UNIX systems, needs testing
    char* dir = const_cast<char*>(filename.c_str());    //convert string to char*
    if (_mkdir(dir) != EEXIST) {
        _mkdir(dir);
    }

    ofile.open(current_path + "//" + filename + "//" + file);
    //Write to file here

    ofile << "Numeric:,";

    if (rel_err != 0) {
        ofile << "Analytic:,";
        ofile << "Relative Error:" << std::endl;
        for (int i = 0; i < n; i++) {
            ofile << v[i] << ",";
            ofile << u[i] << ",";
            ofile << rel_err[i] << std::endl;
        }
    }
    else {
        ofile << "Analytic:" << std::endl;
        for (int i = 0; i < n; i++) {
            ofile << v[i] << ",";
            ofile << u[i] << std::endl;
        }
    }

    ofile.close();
}

double* analyticalSolution(int n, double h) {
    double* u = new double[n];
    for (int i = 0; i < n; i++) {
        u[i] = 1 - (1 - std::exp(-10)) * (i + 1) * h - std::exp(-10 * (i + 1) * h);
    }
    return u;
}

//Begin main program
int main(int argc, char* argv[]) {
    // A*v = h^2 *f_i = b_tilde
    //A needs a, b and c values. 
    //if constant, they could be hard coded
    // v is what we will calculate
    // h is found from n, which are our integration points
    // b_tilde is found from f, which we will calculate

    //Assumed as integers in this case, but we will treat them sort of like arrays
    int a;
    int b;
    int c;

    char task;

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

    askForInput:
        std::cout << "Please input what task you want to solve (valid tasks: b, c, d, e, f): ";
        std::cin >> task;

        std::clock_t start, finish;


        switch (task) {
        case 'b': {
            //TODO: Switch to m to n;
            int n[3] = { 10, 100, 1000 };
            for (int i = 0; i < 3; i++) {
                double h = 1.0 / (n[i] + 1);

                double* v;

                //Time the execution of the algorithm
                start = std::clock();
                v = generalSolver(n[i], h, a, b, c);
                finish = std::clock();

                double* u;
                u = analyticalSolution(n[i], h);

                double execution_time = double(finish - start) / double(CLOCKS_PER_SEC);
                std::cout << "Execution time for n = : " << n[i] << " is: " << std::fixed << execution_time
                    << std::setprecision(7) << " seconds" << std::endl;

                writeToFile("task_b", n[i], v, u);

                delete[] v, u;
                v, u = NULL;
            }
            std::cout << "Task b completed" << std::endl;
            break;
        }
        case 'c': {
            int n[5] = { 10, 100, 1000, 10000, 1000000 };
            std::cout << "Not yet implemented" << std::endl;
            exit(1);
        }
        case 'd': {
            int n = int(pow(10, 7));
            double h = 1.0 / (n + 1);

            double* v;
            double* u;

            //Can use specialized algorithm instead
            v = generalSolver(n, h, a, b, c);
            u = analyticalSolution(n, h);

            //The realtive error
            double* rel_err = new double[n];
            for (int i = 0; i < n; i++) {
                u[i] = 1 - (1 - std::exp(-10)) * (i + 1) * h - std::exp(-10 * (i + 1) * h);
                rel_err[i] = std::log10(std::abs((v[i] - u[i]) / u[i]));
            }

            writeToFile("task_c", n, v, u, rel_err);

            std::cout << "Not completely implemented yet" << std::endl;

            delete[] v, u;
            v, u = NULL;
            break;
        }
        case 'e': {
            std::cout << "Not yet implemented" << std::endl;
            exit(1);
        }
        case 'f': {
            std::cout << "Not yet implemented" << std::endl;
            exit(1);
        }
        case '0': {
            exit(1);
        }
        default: {
            std::cout << "Error, must insert valid character!" << std::endl;
            //Jumps back to cin
            goto askForInput;
        }
        }

    }

    return 0;
}