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

#include "network.h"
#include "coordinate.h"
#include "vessel.h"

#include "Triplet.h"
#include "PressureListCreator.h"
#include "EpetraTutorials.h"

typedef std::vector<Vessel> VesselVector;
typedef std::vector<std::vector<double>> Vector2D;
typedef std::vector<Triplet> TripletVector;

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
	
	Network network(numLevels, dims);
	if (!mpiRank) { 
		std::cout << "Generations: " << numLevels << std::endl;
		std::cout << "Total Number of Vessels: " << (pow(2, numLevels+1) - 2) << std::endl; 
	}

	// Generate network and return vector containing vessels
	double time = MPI_Wtime();
	VesselVector vessels = network.generate(filename);
	if (!mpiRank) {
		std::cout << "-- Total Network Generation Time: " << (MPI_Wtime()-time) << " seconds" << std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// TripletVector TripletList;
	// CreatePressureTripletList(filename, numLevels, TripletList);

	// Epetra tests
	//MLAztecOO ();
	if (mpiSize == 1) {
		AssembleMatrixSequential(network, filename, numLevels);
	}
	

	MPI_Finalize();
	return 0;
}	