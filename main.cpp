#include <iostream>
#include "time.h"
#include "main.h"

int main()
{
  /*
  Program for solving the linear equation Av = b_tilde, when A is a tridiagonal matrix.
  This program assumes that it will be given a single value each for a,b, and c, the sub-, super-diagonals and the diagonal, however, 
  it would be trivial to add the functionality to read in an array of values.
  */

  int a = -1;
  int b = 2;
  int c = -1;
  char task = ' ';

  std::clock_t start, finish;

  //This part uses a while loop in conjunction with a switch to keep the program running until
  //the user wants to exit
  std::cout << "Valid tasks are b, c, d, e, f and 0 to exit. " << std::endl;
  while (task != '0')
  {
    std::cout << "Insert task: ";
    std::cin >> task;

    switch (task)
    {
    case 'b':
    {
      int n[3] = {10, 100, 1000};
      for (int i = 0; i < 3; i++)
      {
        double h = 1.0 / (n[i] + 1);

        //Timing the execution of the algorithm
        start = std::clock();
        double *v = generalSolver(n[i], h, a, b, c);
        finish = std::clock();

        double *u = analyticalSolution(n[i], h);

        double execution_time = double(finish - start) / double(CLOCKS_PER_SEC);
        std::cout << "Execution time for n = " << n[i] << " is " << std::fixed << std::setprecision(4)
                  << execution_time * 1000 << "ms" << std::endl;

        writeToFile("task_b", n[i], v, u);
        writeExecTimeToFile("task_b", n[i], execution_time);

        delete[] v, u;
        v, u = NULL;
      }
      std::cout << "Task b has been completed \n"
                << std::endl;
      break;
    }
    case 'c':
    {
      int n[6] = {10, 100, 1000, 10000, 100000, 1000000};
      for (int i = 0; i < 6; i++)
      {

        double h = 1.0 / (n[i] + 1);

        start = std::clock();
        double *v_spec = specSolver(n[i], h);
        finish = std::clock();

        double execution_time_spec = double(finish - start) / double(CLOCKS_PER_SEC);
        std::cout << "Execution time for the special case, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                  << execution_time_spec * 1000 << "ms" << std::endl;

        //comparing with the general solution
        //includes writing the analytical expression for b_v in the special case in the timing, which the lecturer said shouldn't really count
        start = std::clock();
        double *v_gen = generalSolver(n[i], h, a, b, c);
        finish = std::clock();

        double execution_time_gen = double(finish - start) / double(CLOCKS_PER_SEC);
        std::cout << "Execution time for the general case, for n = : " << n[i] << " is: " << std::fixed << std::setprecision(4)
                  << execution_time_gen * 1000 << "ms" << std::endl;

        double *u = analyticalSolution(n[i], h);

        writeToFile("task_c", n[i], v_spec, u);
        writeExecTimeToFile("task_c", n[i], execution_time_spec);

        delete[] v_spec, v_gen, u;
        v_spec, v_gen, u = NULL;
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
        v, u, rel_err = NULL;
      }
      delete[] rel_err_max;
      rel_err_max = NULL;

      break;
    }
    case 'e':
    {
      int n[3] = {10, 100, 1000}; //10 000 works, about 33 seconds, but 100 000 runs out of memory
      int m = 50;                 //number of executions to average over
      for (int i = 0; i < 3; i++)
      {

        double h = 1.0 / (n[i] + 1);
        double *execution_time = new double[m];
        double execution_time_avg = 0; //finding average execution time over 10 tries

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
        v_lu, u = NULL;
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
      double execution_time_avg_lu = 0;
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

  return 0;
}
