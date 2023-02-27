
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

	long long nnz = rowArray[numVessels];
    output.adjncy = (idx_t*) malloc(nnz * sizeof(idx_t));
    output.nnz = nnz;
    idx_t* colArray = output.adjncy;

	for (size_t i = 0; i < numVessels; i++)
    {
        const std::vector<long long>& connectedVessels = vessels[i].getConnectedVessels();
        idx_t* colCopyPosition = colArray + rowArray[i];
        std::copy(connectedVessels.begin(), connectedVessels.end(), colCopyPosition);
    }

    real_t uniformWeight = 1.0 / static_cast<double>(numParts);
    std::fill(output.tpwgts, output.tpwgts + numParts, uniformWeight);

    std::vector<long long> numVesselsPerNode(mpiSize);

	MPI_Allgather(&numVessels, 1, MPI_LONG_LONG_INT, numVesselsPerNode.data(), 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
    output.vtxdist[0] = 0;
    std::partial_sum(numVesselsPerNode.begin(), numVesselsPerNode.end(), output.vtxdist + 1);

    *output.ubvec = 1.08;
}