
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "PressureListCreator.h"
#include "Triplet.h"
#include "VesselGenerator/vessel.h"
#include "VesselGenerator/ImportExportCommon.h"
#include "VesselGenerator/Hdf5FileStructureNames.h"
#include "VesselGenerator/ImportExportHdfBlock.h"

typedef std::vector<Vessel> VesselVector;
typedef std::vector<Triplet> TripletVector;

// void ImportNodesAndConductances(std::string filename, long long numVessels, std::vector<long long>& nodeArray, std::vector<double>& conductanceArray) {
// 	hid_t fileId = OpenHdfFile(filename);
	
// 	hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
// 	CollectiveHdfBlockImport<long long, 2, 2>(nodeDatasetId, 0, numVessels, nodeArray);
// 	H5Dclose(nodeDatasetId);

// 	hid_t conductanceDatasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, HEALTHY_CONDUCTANCE_DATASET);
// 	CollectiveHdfBlockImport<double, 2, 1>(conductanceDatasetId, 0, numVessels, conductanceArray);
// 	H5Dclose(conductanceDatasetId);

// 	H5Fclose(fileId);
// }

void ImportNodesAndConductances(std::string filename, int start, int numRows, std::vector<long long>& nodeArray, std::vector<double>& conductanceArray) {
	hid_t fileId = OpenHdfFile(filename);
	
	hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
	CollectiveHdfBlockImport<long long, 2, 2>(nodeDatasetId, start, numRows, nodeArray);
	H5Dclose(nodeDatasetId);

	hid_t conductanceDatasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, HEALTHY_CONDUCTANCE_DATASET);
	CollectiveHdfBlockImport<double, 2, 1>(conductanceDatasetId, start, numRows, conductanceArray);
	H5Dclose(conductanceDatasetId);

	H5Fclose(fileId);
}

void CreatePressureTripletList(std::string filename, int numLevels, TripletVector& TripletList) {
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	long long numVessels = pow(2, numLevels+1) - 2;

	// determine number of rows per node
	int count = numVessels / mpiSize;
	int remainder = numVessels % mpiSize;
	int start = mpiRank * count + std::min(mpiRank, remainder);
	int end = (mpiRank + 1) * count + std::min(mpiRank + 1, remainder);
	int numRows = end-start;

	//std::cout << "Rank: " << mpiRank << ", start = " << start << ", numRows = " << numRows << std::endl;
	
	// Import Arrays
	std::vector<long long> nodeArray;
	std::vector<double> conductanceArray;
	ImportNodesAndConductances(filename, start, numRows, nodeArray, conductanceArray);

	//! print conductances and nodes
	// std::cout << "Values on Rank: " << mpiRank << std::endl;
	// int index = 0;
	// for(int i = 0; i < conductanceArray.size(); i++) {
	// 	std::cout << conductanceArray[i] << std::endl;
	// 	std::cout << nodeArray[index] << ", " << nodeArray[index+1] << std::endl;
	// 	index += 2;
	// }

	//! assemble triplet list
	int index=0;
	int mat_count = 0;
	for (int i = 0; i < numRows; ++i) {
		double conductance = conductanceArray[i];
		const int node1 = nodeArray[index];
		const int node2 = nodeArray[index+1];

		if (node1 == 0) {
			TripletList.push_back(Triplet(node1, node1, 1.0f));
			mat_count +=1;
		} else {
			TripletList.push_back(Triplet(node1, node1, -conductance));
			TripletList.push_back(Triplet(node1, node2, conductance));
			mat_count +=2;
		}
		if (node2 == 0) {
			TripletList.push_back(Triplet(node2, node2, 1.0f));
			mat_count +=1;
		} else {
			TripletList.push_back(Triplet(node2, node2, -conductance));
			TripletList.push_back(Triplet(node2, node1, conductance));
			mat_count +=2;
		}
		index += 2;
	}

	//std::cout << "Count = " << mat_count << ", on rank " << mpiRank << std::endl;

	//!print triplet list
	// for (int i = 0; i < TripletList.size(); i++) {
	// 	std::cout << TripletList[i] << std::endl;
	// }
}