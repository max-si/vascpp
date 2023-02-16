
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include <mpi.h>
#include <parmetis.h>

#include "BloodFlow/Triplet.h"
#include "AdjacencyMatrix.h"
#include "PartitionFunctions.h"
#include "PreprocessorVessel.h"

//typedef PreprocessorVessel Vessel;

//! This adjacency matrix maps the graph edges, not vertices*
void GenerateAdjacencyMatrix(std::vector<PreprocessorVessel>& vessels, int numParts, AdjacencyMatrix& output) {
	int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	// if(!mpiRank) {
	// 	for(int i=0; i < vessels.size(); i++) {
	// 		std::cout << vessels[i] << std::endl;
	// 	}
	// }

	size_t numVessels = vessels.size();

	output.AllocateMatrix(numParts, numVessels);

	idx_t* rowArray = output.xadjx;

	rowArray[0] = 0;

	for (size_t i = 0; i < numVessels; i++)
    {
        rowArray[i + 1] = vessels[i].getNumConnectedVessels();
    }

	for (size_t i = 1; i < numVessels + 1; i++)
    {
        rowArray[i] += rowArray[i - 1];
    }

	std::cout << mpiRank << " Row Array" << std::endl;
	for(int i = 0; i < numVessels+1; i++) {
		std::cout << rowArray[i] << ", ";
	} std::cout << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);

	long long nnz = rowArray[numVessels];
	std::cout << "nnz = " << nnz << ", mpiRank = " << mpiRank << std::endl;
    output.adjncy = (idx_t*) malloc(nnz * sizeof(idx_t));
    output.nnz = nnz;
    idx_t* colArray = output.adjncy;

	for (size_t i = 0; i < numVessels; i++)
    {
        const std::vector<long long>& connectedVessels = vessels[i].getConnectedVessels();
        idx_t* colCopyPosition = colArray + rowArray[i];
        std::copy(connectedVessels.begin(), connectedVessels.end(), colCopyPosition);
    }

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << mpiRank << " Col Array" << std::endl;
	for(int i = 0; i < nnz; i++) {
		std::cout << colArray[i] << ", ";
	} std::cout << std::endl;


	//! TRY 2D FIRST
	// size = d*n_i where n_i is number of locally stored vertices and d is dimensions
	//real_t* vertxCoords = output.xyz;
	//output.xyz = (real_t*) malloc(nnz * sizeof(idx_t));
	// real_t x, y;
	// if(!mpiRank)
	// 	for (size_t i = 0; i < numVessels; i++) {
	// 		x = vessels[i].get_start_x();
	// 		y = vessels[i].get_start_y();
	// 		std::cout << x << ", " << y << std::endl;
	// 	}


    real_t uniformWeight = 1.0 / static_cast<double>(numParts);
    std::fill(output.tpwgts, output.tpwgts + numParts, uniformWeight);

    std::vector<long long> numVesselsPerNode(mpiSize);

	MPI_Allgather(&numVessels, 1, MPI_LONG_LONG_INT, numVesselsPerNode.data(), 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
    output.vtxdist[0] = 0;
    std::partial_sum(numVesselsPerNode.begin(), numVesselsPerNode.end(), output.vtxdist + 1);

	// if(!mpiRank)
	// 	for (int i = 0; i < mpiSize; i++) {
	// 		std::cout << numVesselsPerNode[i] << std::endl;
	// 	}

    *output.ubvec = 1.08;

	//std::cout << << std::endl;
}