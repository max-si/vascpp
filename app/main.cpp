// vascpp main file
// TODO: clean up files/folders, add comments, etc.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <mpi.h>

#include "VesselGenerator/network.h"
#include "VesselGenerator/coordinate.h"
#include "VesselGenerator/vessel.h"

#include "NetworkPreprocessing/PartitionFunctions.h"

#include "BloodFlow/BloodFlowFunc.h"
#include "BloodFlow/SequentialSolver.h"

typedef std::vector<Vessel> VesselVector;
typedef std::vector<std::vector<double>> Vector2D;

int main(int argc, char* argv[]) {
	// Initialize network(levels, dimensions) from user input and define filename
	std::string filename;
	int numLevels, dims;
    if (argc < 2)
    {
		filename = "test";
        numLevels = 12;
        dims = 2;
    }
    else if (argc == 2)
    {
		filename = "test";
        numLevels = atoi(argv[1]);
        dims = 2;
    }
    else if (argc == 3)
    {
		filename = "test";
        numLevels = atoi(argv[1]);
        dims = atoi(argv[2]);
    }
    else if (argc == 4)
    {
        numLevels = atoi(argv[1]);
        dims = atoi(argv[2]);
		filename = argv[3];
    }
    else
    {
        std::cout << "incorrect input arguments" << std::endl;
		MPI_Finalize();
        return 1;
    }	
	filename += ".h5";
	
	MPI_Init(NULL, NULL);
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	if (!mpiRank) { 
		std::cout << R"(__        __   ____     ______     ______   _______    ________)" << std::endl;
		std::cout << R"(\+\      /+/  /+/\+\   /+/`````   /+/````` |+|````\+\ |+|````\+\)" << std::endl;
		std::cout << R"( \+\    /+/  /+/  \+\  \+\_____  /+/       |+|,,,,/+/ |+|,,,,/+/)" << std::endl;
		std::cout << R"(  \+\  /+/  /+/____\+\       \+\ \+\       |+|        |+|)" << std::endl;
		std::cout << R"(   \+\/+/  /+/      \+\ ,,,,,/+/  \+\,,,,, |+|        |+|)" << std::endl;
		std::cout << "    \u203E\u203E\u203E\u203E   \u203E\u203E        \u203E\u203E \u203E\u203E\u203E\u203E\u203E\u203E\u203E";
		std::cout << "    \u203E\u203E\u203E\u203E\u203E\u203E\u203E \u203E\u203E\u203E        \u203E\u203E\u203E" << std::endl;
	}
	
	Network network(numLevels, dims);
	if (!mpiRank) { 
		std::cout << std::endl << "-----------------------------" << std::endl;
		std::cout << "- Building Network Geometry -" << std::endl;
		std::cout << "-----------------------------" << std::endl;
		std::cout << "Generations: " << numLevels << std::endl;
		std::cout << "Total Number of Vessels: " << (pow(2, numLevels+1) - 2) << std::endl; 
	}

	// Generate network and return vector containing vessels
	double time = MPI_Wtime();
	network.generate(filename);
	if (!mpiRank) {
		std::cout << "-- Total Network Generation Time: " << (MPI_Wtime()-time) << " seconds" << std::endl;
	}

	// single node Blood Flow Solver
	if (mpiSize == 1) {
		double solve_time = MPI_Wtime();
		AssembleMatrixSequential(network, filename, numLevels);
		std::cout << "-- Total Sequential Blood Flow Solve Time: " << (MPI_Wtime()-solve_time) << " seconds" << std::endl;
	}


	//************************\\
	//*    Parallel Solver   *\\
	//************************\\

	// Starting PreProcessing
	if (mpiSize != 1) { if(!mpiRank) {std::cout << std::endl << "----------------------------------" << std::endl; }}
	if (mpiSize != 1) { if(!mpiRank) {std::cout << "- Starting Network Preprocessing -" << std::endl; }}
	if (mpiSize != 1) { if(!mpiRank) {std::cout << "----------------------------------" << std::endl; }}
	double t2 = MPI_Wtime();
	if (mpiSize != 1) { networkPreproc(filename); }
	if (mpiSize != 1) { if (!mpiRank) { std::cout << std::endl << "-- Total Network Preprocessing Time: " << (MPI_Wtime()-t2) << " seconds" << std::endl; }}

	// Starting Blood Flow Solver (parallel)
	if (mpiSize != 1) { if(!mpiRank) { std::cout << std::endl << "------------------------------" << std::endl; }}
	if (mpiSize != 1) { if(!mpiRank) { std::cout << "- Starting Blood Flow Solver -" << std::endl; } }
	if (mpiSize != 1) { if(!mpiRank) { std::cout << "------------------------------" << std::endl; }}
	double t3 = MPI_Wtime();
	if (mpiSize != 1) { BloodFlowExe(filename); }
	if (mpiSize != 1) { if (!mpiRank) { std::cout << std::endl << "-- Total Blood Flow App Time: " << (MPI_Wtime()-t3) << " seconds" << std::endl; }}

	MPI_Finalize();
	return 0;
}	