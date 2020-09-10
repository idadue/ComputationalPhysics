#include <iostream>
#include "main.h"

int main(int argc, char *argv[])
{
  /*
  Program for solving the linear equation Av = b_tilde, when A is a tridiagonal matrix.
  This program assumes that it will be given a single value each for a,b, and c, the sub-, super-diagonals and the diagonal, however, 
  it would be trivial to add the functionality to read in an array of values.
  */

  int a;
  int b;
  int c;

  if (argc < 4)
  {
    //If no commands are specified, exit program
    std::cout << "Bad Usage: " << argv[0] << " is missing arguments" << std::endl;
    exit(1);
  }
  else
  {
    a = std::atoi(argv[1]);
    b = std::atoi(argv[2]);
    c = std::atoi(argv[3]);
    char task = ' ';
    std::clock_t start, finish;

    //This part uses a while loop in conjunction with a switch to keep the program running until
    //the user wants to exit
    std::cout << "Valid tasks are b, c, d, e, f and 0 to exit. " << std::endl;
    while (task != '0')
    {
      std::cout << "Insert task: ";
      std::cin >> task;
      int n[6] = {10, 100, 1000, 10000, 100000, 1000000};

      switch (task)
      {
      case 'b':
      {
        for (int i = 0; i < 3; i++)
        {
          time_and_write(generalSolver, n[i], a, b, c, "task_b");
        }
        break;
      }
      case 'c':
      {
        for (int i = 0; i < 6; i++)
        {
          //special algorithm
          time_and_write(generalSolver, n[i], 0, 0, 0, "task_c_spec");
          //general algorithm
          time_and_write(generalSolver, n[i], a, b, c, "task_c_general");
        }
        break;
      }
      case 'd':
      {
        double *rel_err_max = new double[7];

        for (int j = 0; j < 7; j++)
        {
          int n = int(pow(10, j + 1));
          double h = 1.0 / (n + 1);

          double *v = specSolver(n, h);
          double *u = analyticalSolution(n, h);

          //The relative error
          double *rel_err = new double[n];
          rel_err_max[j] = -10000;
          for (int i = 0; i < n; i++)
          {
            rel_err[i] = std::log10(std::abs((v[i] - u[i]) / u[i]));
            if (rel_err_max[j] < rel_err[i])
            {
              rel_err_max[j] = rel_err[i];
            }
          }

          std::cout << "Maximum relative error for n = " << n << " is "
                    << std::fixed << std::setprecision(4) << rel_err_max[j] << std::endl;

          delete[] v, u, rel_err;
        }
        delete[] rel_err_max;

        break;
      }
      case 'e':
      {
        int m = 50; //number of executions to average over
        for (int i = 0; i < 3; i++)
        {

          double h = 1.0 / (n[i] + 1);
          double *execution_time = new double[m];
          double execution_time_avg = 0;

          for (int j = 0; j < m; j++)
          {

            start = std::clock();
            double *v_lu = lusolver(n[i], h, a, b, c);
            finish = std::clock();

            execution_time[j] = double(finish - start) / double(CLOCKS_PER_SEC);
            std::cout << "Execution time using LU decomposition, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                      << execution_time[j] * 1000 << "ms" << std::endl;

            execution_time_avg += (execution_time[j]) / m;
          }

          std::cout << "Average Execution time using LU decomposition, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                    << execution_time_avg * 1000 << "ms" << std::endl;

          double *v_lu = lusolver(n[i], h, a, b, c);
          double *u = analyticalSolution(n[i], h);

          writeToFile("task_e", n[i], v_lu, u);
          writeExecTimeToFile("exec_time_lu", n[i], execution_time_avg);

          delete[] v_lu, u;
        }
        break;
      }
      case 'f':
      {
        int n;
        int m; //number of executions to average over
        std::cout << "Enter value for n (integration points), m (number of executions): ";
        std::cin >> n >> m;
        std::cout << "\n";
        std::cout << "You entered: " << n << ", " << m << " number of executions" << std::endl;

        double h = 1.0 / (n + 1);
        double *execution_time = new double[m];
        double execution_time_avg_spec = 0;
        double execution_time_avg_gen = 0;

        //Timing the execution of each algorithm
        //general case:
        for (int i = 0; i < m; i++)
        {

          start = std::clock();
          generalSolver(n, h, a, b, c);
          finish = std::clock();

          execution_time[i] = double(finish - start) / double(CLOCKS_PER_SEC);
          std::cout << "Execution time, general case, n = " << n << " is " << std::fixed << std::setprecision(4)
                    << execution_time[i] * 1000 << "ms" << std::endl;

          execution_time_avg_gen += execution_time[i] / m;
        }
        std::cout << "Average Execution time, general case, n = " << n << " is " << std::fixed << std::setprecision(4)
                  << execution_time_avg_gen * 1000 << "ms" << std::endl;

        //special case:
        for (int i = 0; i < m; i++)
        {

          start = std::clock();
          specSolver(n, h);
          finish = std::clock();

          execution_time[i] = double(finish - start) / double(CLOCKS_PER_SEC);
          std::cout << "Execution time, special case, n =" << n << " is " << std::fixed << std::setprecision(4)
                    << execution_time[i] * 1000 << "ms" << std::endl;

          execution_time_avg_spec += execution_time[i] / m;
        }
        std::cout << "Average Execution time, special case, n = " << n << " is " << std::fixed << std::setprecision(4)
                  << execution_time_avg_spec * 1000 << "ms" << std::endl;

        writeExecTimeToFile("exec_time_spec", n, execution_time_avg_spec);
        writeExecTimeToFile("exec_time_gen", n, execution_time_avg_gen);

        std::cout << "Task completed! \n"
                  << std::endl;

        delete[] execution_time;

        break;
      }
      default:
      {
        if (task == '0')
        {
          std::cout << "Exiting..." << std::endl;
        }
        else
        {
          std::cout << "Error, must insert valid character!" << std::endl;
          break;
        }
      }
      }
    }
  }
  return 0;
}
