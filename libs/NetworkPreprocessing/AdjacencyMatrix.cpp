
#include <cstdlib>

#include <mpi.h>
#include <parmetis.h>

#include "AdjacencyMatrix.h"

AdjacencyMatrix::AdjacencyMatrix(idx_t npartsIn, idx_t numVessels, idx_t numEntries)
{
	AllocateMatrix(npartsIn, numVessels, numEntries);
}

AdjacencyMatrix::~AdjacencyMatrix()
{
	DeallocateMatrix();
}

void AdjacencyMatrix::AllocateMatrix(idx_t npartsIn, idx_t numVessels, idx_t numEntries)
{
	if(isAllocated) {
		DeallocateMatrix();
		isAllocated = false;
	}
	nparts = npartsIn;
	numRows = numVessels;
	nnz = numEntries;
	int numProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	vtxdist = (idx_t*) malloc((numProcs+1) * sizeof(idx_t));
	xadjx = (idx_t*) malloc((numVessels + 1) * sizeof(idx_t));
    adjncy = (idx_t*) malloc(numEntries * sizeof(idx_t));
    options = (idx_t*) malloc(4 * sizeof(idx_t));
    options[0] = 1;
    options[1] = 0;
    options[2] = 133829837;
    options[3] = PARMETIS_PSR_COUPLED;
    partArray = (idx_t*) malloc(numVessels * sizeof(idx_t));
    tpwgts = (real_t*) malloc(nparts * sizeof(real_t));
    ubvec = (real_t*) malloc(sizeof(real_t));
    isAllocated = true;
}

void AdjacencyMatrix::AllocateMatrix(idx_t npartsIn, idx_t numVessels)
{
    if (isAllocated)
    {
        DeallocateMatrix();
        isAllocated = false;
    }
    nparts = npartsIn;
    numRows = numVessels;
    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    vtxdist = (idx_t*) malloc((numProcs + 1) * sizeof(idx_t));
    xadjx = (idx_t*) malloc((numVessels + 1) * sizeof(idx_t));
    options = (idx_t*) malloc(4 * sizeof(idx_t));
    options[0] = 1;
    options[1] = 0;
    options[2] = 133829837;
    options[3] = PARMETIS_PSR_COUPLED;
    partArray = (idx_t*) malloc(numVessels * sizeof(idx_t));
    tpwgts = (real_t*) malloc(nparts * sizeof(real_t));
    ubvec = (real_t*) malloc(sizeof(real_t));
    isAllocated = true;
}

void AdjacencyMatrix::DeallocateMatrix()
{
    if (isAllocated)
    {
        free(vtxdist);
        free(xadjx);
        free(adjncy);
        free(options);
        free(partArray);
        free(tpwgts);
        free(ubvec);
        free(vwgt);
        free(adjwgt);
        isAllocated = false;
    }
}

AdjacencyMatrix::AdjacencyMatrix(){};