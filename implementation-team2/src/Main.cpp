/*
 * Main.cpp
 *
 *  Created on: May 7, 2015
 *      Authors: Sarah Lutteropp, Pierre Barbera, Marin Vlastelica Pogančić
 */

#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <algorithm>
#include "TKF91Sequential.h"

//Default parameter values
const double DEFAULT_LAMBDA = 1, DEFAULT_MU = 2;
const double DEFAULT_TAU = 0.1;
const double DEFAULT_PI[4] = {0.25, 0.25, 0.25, 0.25};
long double benchmark_method(std::string&, std::string&);

long double benchmark_method(std::string& seq1, std::string& seq2)
{
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 10; ++i)
            TKF91Sequential(seq1, seq2, DEFAULT_LAMBDA, DEFAULT_MU, DEFAULT_TAU, DEFAULT_PI);
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long double ms =  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        ms /= 10.;
        std::cout << "Seq1 length: " << seq1.length() <<  "  " << "Seq2 length: " << seq2.length() << std::endl;
        std::cout << "Elapsed Time in millisecs: " << ms << std::endl;
    	return ms;
}


int main(int argc, char** argv)
{
    //Default parameter values
    std::string inputA = "", inputB = "";
    double inputLambda = DEFAULT_LAMBDA, inputMu = DEFAULT_MU;
    double inputTime = DEFAULT_TAU;
    double inputPi[4] = {0.25, 0.25, 0.25, 0.25};
    bool outputMatrices = false;
    bool outputAlignment = false;
    bool benchmark = false;
    int runs = 1;

    for(int i = 0; i < argc; i++)
        if (strcmp(argv[i], "--output-matrices") == 0) {
            outputMatrices = true;
        } else if (strcmp(argv[i], "--output-alignment") == 0) {
            outputAlignment = true;
        } else if (strcmp(argv[i], "-t") == 0) {
            benchmark = true;
        } else if (strcmp(argv[i], "-h") == 0) {
            argc = 0;
        } 


    if(argc != 0)
    {
        for (int i = 1; i < argc; i += 2) {
            if (i + 1 != argc) {
                if (strcmp(argv[i], "--sequence-1") == 0) {
                    inputA = argv[i + 1];
                } else if (strcmp(argv[i], "--sequence-2") == 0) {
                    inputB = argv[i + 1];
                } else if (strcmp(argv[i], "--lambda") == 0) {
                    inputLambda = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--mu") == 0) {
                    inputMu = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--tau") == 0) {
                    inputTime = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--pa") == 0) {
                    inputPi[0] = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--pc") == 0) {
                    inputPi[1] = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--pg") == 0) {
                    inputPi[2] = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "--pt") == 0) {
                    inputPi[3] = atof(argv[i + 1]);
                } else if (strcmp(argv[i], "-r") == 0) {
                    runs = atoi(argv[i + 1]);
                }
            }
        }
    }
    if(argc == 0)
    {

        //Print default values
        std::cout << "Default lambda value: \t" << inputLambda << std::endl;
        std::cout << "Default mu value: \t" << inputMu << std::endl;
        std::cout << "Default tau value: \t" << inputTime << std::endl;
        std::cout << "Default nucleotide probabilities: \t" << inputPi[0] << " " << inputPi[1] << " "
                  << inputPi[2] << " " <<  inputPi[3] << " " << "\n\n";

        //Usage explanation
        std::cout << "-----------USAGE------------" << "\n\n";
        std::cout << "<prog> --sequence-1 <SEQUENCE-1> --sequence-2 <SEQUENCE-2> [--lambda <LAMBDA>]"
                     " [--mu <MU>] [--tau <TAU>] [--pa <PA>] [--pc <PC>] [--pg <PG>] [--pt <PT>] [-t] [-v] [-r] [-h]" << "\n\n";
        std::cout << "SEQUENCE-1\t The first nucleotide sequence." << std::endl;
        std::cout << "SEQUENCE-2\t The second nucleotide sequence." << std::endl;
        std::cout << "LAMBDA\t\t The birth rate." << std::endl;
        std::cout << "MU\t\t The death rate." << std::endl;
        std::cout << "TAU\t\t The time." << std::endl;
        std::cout << "PA\t\t Probability for nucleotide A." << std::endl;
        std::cout << "PC\t\t Probability for nucleotide C." << std::endl;
        std::cout << "PG\t\t Probability for nucleotide G." << std::endl;
        std::cout << "PT\t\t Probability for nucleotide T." << std::endl;
        std::cout << "--output-matrices\t\t Print matrices." << std::endl;
        std::cout << "--output-alignment\t\t Print alignment." << std::endl;
        std::cout << "-t\t\t Benchmark data." << std::endl;
        std::cout << "-r\t\t Specify number of times kernel execution is repeated." << std::endl;
        exit(0);

    }
    else if((inputA == "" || inputB == "") && !benchmark)
    {
        std::cout << "Incorrect usage, use flag -h for help!" << std::endl;
        exit(1);
    }



    if(benchmark)
    {
        std::ifstream input("../doc/sequences");
        std::ofstream output;
        output.open("../doc/benchmark.csv", std::ofstream::app);
        if(!output.is_open())
        {
            std::cerr << "Failed to open output file. " << std::endl;
            exit(-1);
        }
        if(input.is_open())
        {
            std::string seq1, seq2;
            size_t i = 0;
            long double times[4];
            while(std::getline(input, seq1))
            {
                std::getline(input, seq2);
            
                std::replace(seq1.begin(), seq1.end(), 'N', 'A');
                std::replace(seq2.begin(), seq2.end(), 'N', 'A');

                std::cout << "Doing benchmark " << i << std::endl;
                times[i] = benchmark_method(seq1, seq2);
                i+=1;

            }
            output << "log_caching_round_up," << times[0] << "," << times[1] << "," << times[2]
            << "," << times[3] <<  std::endl;
            output.close();
            input.close();
            std::cout << "Benchmarked succesfully!" << std::endl;
            exit(0);
        }
        else
        {
            std::cerr<< "Failed to open input file." << std::endl;
            exit(-1);
        }
    }
   	else
   	{
   		std::cout << "log_caching_round_up"<< std::endl;
        std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < runs; ++i)
    	   TKF91Sequential(inputA, inputB, inputLambda, inputMu, inputTime, inputPi, outputMatrices, outputAlignment);
        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        long double ms =  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        ms /= double(runs);
        std::cout << "Average time (ms): " << ms << std::endl;
        
	}
    return 0;
}

