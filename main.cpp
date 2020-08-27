#include <iostream>
#include <fstream>
#include <string>

//include setw, setprecision, etc
#include <iomanip>


//begin main program
int main(int argc, char *argv[]){

    //Input from command line
    // A*u = h^2 b_tilde
    //A needs a, b and c values. Could be vector of constants, but need only be single constant for this project?
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
    
    //Append filename with value of n, for example filename_n.csv
    file.append("_" + std::to_string(n));

    //Calculation starts here
     
    //Calculation ends here

    //Open and close file
    ofile.open(file);
    //Write to file here
    //can format the output data later

    ofile.close();


    return 0;
}