#include <iostream>
#include "time.h"
#include "main.h"

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
                writeToFile("task_d", n, v, u, rel_err);

                std::cout << "Writing completed \n" << std::endl;
                //std::cout << "Not completely implemented yet" << std::endl;

                delete[] v, u, rel_err;
                v, u, rel_err = NULL;
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