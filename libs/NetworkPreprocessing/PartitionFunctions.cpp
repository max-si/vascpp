
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include <mpi.h>
#include <parmetis.h>

#include "VesselGenerator/ImportExportCommon.h"
#include "preImporter.h"
#include "AdjacencyMatrix.h"
#include "PartitionFunctions.h"
#include "SortedVesselExporter.h"
#include "NodeProcessing.h"


void networkPreproc(std::string filename) {
    //std::cout << "Entering network preprocessor" << std::endl;
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	hid_t fileId = OpenHdfFile(filename);
	std::vector<PreprocessorVessel> vessels = PreprocessorVesselImporter(fileId);
	MPI_Barrier(MPI_COMM_WORLD);

	int numParts = mpiSize;
	AdjacencyMatrix adjacencyMatrix;
	GenerateAdjacencyMatrix(vessels, numParts, adjacencyMatrix);

	PartitionMatrix(vessels, adjacencyMatrix);

    adjacencyMatrix.DeallocateMatrix();

    SortVesselsByPartition(vessels);

    RenumberConnectedVesselsBeforeExport(vessels);

    ExportVesselDataSorted(fileId, vessels, numParts);

    // if(!mpiRank)
    //     for(int i = 0; i < vessels.size(); i++) {
    //         std::cout << vessels[i]  << std::endl;
    //     }

    NodeProcessing(fileId);

    //std::cout << "Closing File" << std::endl;
    H5close();
}


void PartitionMatrix(std::vector<PreprocessorVessel>& vessels, AdjacencyMatrix& adjacencyMatrix)
{
    idx_t* vtxdist = adjacencyMatrix.vtxdist;
    idx_t* xadjx = adjacencyMatrix.xadjx;
    idx_t* adjncy = adjacencyMatrix.adjncy;
    idx_t* nparts = &adjacencyMatrix.nparts;
    idx_t* edgecut = &adjacencyMatrix.edgecut;
    idx_t* options = adjacencyMatrix.options;
    idx_t* partArray = adjacencyMatrix.partArray;

    idx_t* ncon = &adjacencyMatrix.ncon;
    real_t* tpwgts = adjacencyMatrix.tpwgts;
    real_t* ubvec = adjacencyMatrix.ubvec;
    idx_t* numflag = &adjacencyMatrix.numflag;
    idx_t* wgtflag = &adjacencyMatrix.wgtflag;
    idx_t* vwgt = adjacencyMatrix.vwgt;
    idx_t* adjwgt = adjacencyMatrix.adjwgt;

    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //if (mpiSize ==1)
    //TestForSymmetry(adjacencyMatrix);

    //abort();

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    //cout << "entering parMetis" << endl;

    int errorVal = ParMETIS_V3_PartKway(vtxdist, xadjx, adjncy, vwgt, adjwgt,
        wgtflag, numflag, ncon, nparts, tpwgts, ubvec, options, edgecut,
        partArray, &comm);

    if (errorVal != METIS_OK)
    {
        int myRank;
        MPI_Comm_rank(comm, &myRank);
        if (myRank == 0)
            std::cout << "Error Occurred: Running ParMetis with debug output "
                         "for more details"
                      << std::endl;
        options[0] = 1;
        options[1] = PARMETIS_DBGLVL_INFO;
        errorVal = ParMETIS_V3_PartKway(vtxdist, xadjx, adjncy, vwgt, adjwgt,
            wgtflag, numflag, ncon, nparts, tpwgts, ubvec, options, edgecut,
            partArray, &comm);
        exit(1);
    }
    if (!mpiRank)
        std::cout << "Initial Edge Count: " << *edgecut << std::endl;
    long long i, numVessels = vessels.size();
    for (i = 0; i < numVessels; i++)
    {
        vessels[i].setPartitionNumber(partArray[i]);
    }
}

void SortVesselsByPartition(std::vector<PreprocessorVessel>& vessels)
{
    //pss::parallel_sort(vessels.begin(), vessels.end(), CompareVesselNodesPartition);
    std::sort(vessels.begin(), vessels.end(), CompareVesselNodesPartition);
}