#include <iostream>
#include "time.h"
#include "main.h"

//TODO Find also the precise number of floating point operations needed to solve equations.

//Begin main program
int main(int argc, char* argv[]) {
    // A*v = h^2 *f_i = b_tilde
    //A needs a, b and c values.
    //if constant, they could be hard coded
    // v is what we will calculate
    // h is found from n, which are our integration points
    // b_tilde is found from f, which we will calculate

    //Assumed as integers in this case, but we will treat them like arrays in the implementation
    //If we wanted to use n different integers for each diagonal we would have to rewrite some of the code, and it would be better to
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

        //This part uses a while loop in conjunction with a switch to keep the program running until
        //the user wants to exit
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
                    writeExecTimeToFile("task_b", n[i], execution_time);

                    delete[] v, u;
                    v, u = NULL;
                }
                std::cout << "Task b has been completed \n" << std::endl;
                break;
            }
            case 'c': {
                //TODO: Create a specialized algo
                //Similar structure as task b
                int n[6] = { 10, 100, 1000, 10000, 100000, 1000000 };
                for (int i = 0; i < 6; i++) {

                  double h = 1.0 / (n[i] + 1);

                  //Time the execution of the algorithm
                  start = std::clock();
                  double* v_spec = specSolver(n[i], h);
                  finish = std::clock();

                  double execution_time_spec = double(finish - start) / double(CLOCKS_PER_SEC);
                  std::cout << "Execution time for the special case, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                  << execution_time_spec*1000 << "ms" << std::endl;

                  //comparing with the general solution
                  //includes writing the analytical expression for b_v in the special case in the timing, which the lecturer said shouldn't really count
                  start = std::clock();
                  double* v_gen = generalSolver(n[i], h, a, b, c);
                  finish = std::clock();

                  double execution_time_gen = double(finish - start) / double(CLOCKS_PER_SEC);
                  std::cout << "Execution time for the general case, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                  << execution_time_gen*1000 << "ms" << std::endl;

                  double* u = analyticalSolution(n[i], h);

                  writeToFile("task_c", n[i], v_spec, u);
                  writeExecTimeToFile("task_c", n[i], execution_time_spec);

                  delete[] v_spec, v_gen, u;
                  v_spec, v_gen, u = NULL;
                }
                break;
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
                int n[3] = { 10, 100, 1000}; //10 000 works, about 33 seconds, but 100 000 runs out of memory
                for ( int i = 0; i < 3; i++){

                  double h = 1.0 / (n[i] + 1);
                  double* execution_time = new double[10];
                  double execution_time_avg = 0; //finding average execution time over 10 tries

                  for ( int j = 0; j < 10; j++){

                    start = std::clock();
                    double* v_lu = lusolver(n[i], h, a, b, c);
                    finish = std::clock();


                    execution_time[j] = double(finish - start) / double(CLOCKS_PER_SEC);
                    std::cout << "Execution time using LU decomposition, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                    << execution_time[j]*1000 << "ms" << std::endl;

                    execution_time_avg += (execution_time[j])/10;


                  }

                std::cout << "Average Execution time using LU decomposition, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                << execution_time_avg*1000 << "ms" << std::endl;

                double* v_lu = lusolver(n[i], h, a, b, c);
                double* u = analyticalSolution(n[i], h);

                writeToFile("task_e", n[i], v_lu, u);
                writeExecTimeToFile("exec_time_lu", n[i], execution_time_avg);



                delete[] v_lu, u;
                v_lu, u = NULL;
                }
                break;
            }
            case 'f': {
                int n; //std::string filename;
                int m; //number of executions to average over
                std::cout << "Enter value for n, m (number of executions): ";
                std::cin >> n >> m;
                std::cout << "\n";
                std::cout << "You entered: " << n << ", " << m << " number of executions"<< std::endl;

                double h = 1.0 / (n + 1);
                double* execution_time = new double[m];
                double execution_time_avg_lu = 0;
                double execution_time_avg_spec = 0;
                double execution_time_avg_gen = 0;

                //Timing the execution of each algorithm

                //general case:
                for ( int i = 0; i < m; i++){

                  start = std::clock();
                  generalSolver(n, h, a, b, c);
                  finish = std::clock();

                  execution_time[i] = double(finish - start) / double(CLOCKS_PER_SEC);
                  std::cout << "Execution time, general case, n = " << n << " is " << std::fixed << std::setprecision(4)
                  << execution_time[i] * 1000 << "ms" << std::endl;

                  execution_time_avg_gen += execution_time[i]/m;

                }
                std::cout << "Average Execution time, general case, n = " << n << " is " << std::fixed << std::setprecision(4)
                << execution_time_avg_gen * 1000 << "ms" << std::endl;

                //special case:
                for ( int i = 0; i < m; i++){

                  start = std::clock();
                  specSolver(n, h);
                  finish = std::clock();

                  execution_time[i] = double(finish - start) / double(CLOCKS_PER_SEC);
                  std::cout << "Execution time, special case, n =" << n << " is " << std::fixed << std::setprecision(4)
                  << execution_time[i] * 1000 << "ms" << std::endl;

                  execution_time_avg_spec += execution_time[i]/m;

                }
                std::cout << "Average Execution time, special case, n = " << n << " is " << std::fixed << std::setprecision(4)
                << execution_time_avg_spec * 1000 << "ms" << std::endl;

                //lu decomposition, too time consuming compared with gen/spec so moved back to e)
                /*
                for ( int i = 0; i < m; i++){

                  start = std::clock();
                  lusolver(n, h, a, b, c);
                  finish = std::clock();

                  execution_time[i] = double(finish - start) / double(CLOCKS_PER_SEC);
                  std::cout << "Execution time, LU-decomp, n =" << n << " is " << std::fixed << std::setprecision(4)
                  << execution_time[i] * 1000 << "ms" << std::endl;

                  execution_time_avg_lu += execution_time[i]/m;

                }
                std::cout << "Average Execution time, LU-decomp, n = " << n << " is " << std::fixed << std::setprecision(4)
                << execution_time_avg_lu * 1000 << "ms" << std::endl;
                */


                //writeExecTimeToFile("exec_time_lu", n, execution_time_avg_lu);
                writeExecTimeToFile("exec_time_spec", n, execution_time_avg_spec);
                writeExecTimeToFile("exec_time_gen", n, execution_time_avg_gen);

                std::cout << "Task completed! \n" << std::endl;

                delete[] execution_time;

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
