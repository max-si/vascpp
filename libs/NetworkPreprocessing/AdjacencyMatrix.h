#pragma once

#include <parmetis.h>

// Adjacency Matrix class
// Description of the adjacency matrix for graph partitioning algorithm
// Handles all information needed for PARMETIS

class AdjacencyMatrix {
public:
	idx_t* vtxdist = nullptr;
	idx_t* xadjx = nullptr;
	idx_t* adjncy = nullptr;
	idx_t nparts = 0;
	idx_t edgecut = 0;
	idx_t* options = nullptr;
	idx_t* partArray = nullptr;
	idx_t* ndims = nullptr;
	real_t* xyz = nullptr;

	idx_t numRows = 0;
	idx_t nnz = 0;

	idx_t ncon = 1;
	real_t* tpwgts = nullptr;
	real_t* ubvec = nullptr;
	idx_t numflag = 0;
	idx_t wgtflag = 0;
	idx_t* vwgt = nullptr;
	idx_t* adjwgt = nullptr;

	AdjacencyMatrix(idx_t nparts, idx_t numVessels, idx_t numEntries);
	AdjacencyMatrix();
	~AdjacencyMatrix();

	void AllocateMatrix(idx_t npartsIn, idx_t numVessels, idx_t numEntries);
	void AllocateMatrix(idx_t npartsIn, idx_t numVessels);
	void DeallocateMatrix();

private:
	bool isAllocated = false;
};