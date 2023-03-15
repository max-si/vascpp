// This file contains the Network class and methods
// The class generates the vascular network based on instances of Vessel
// TODO: Add comments from pyvascular

#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>

#include <hdf5.h>
#include <mpi.h>

#include "network.h"
#include "coordinate.h"
#include "vessel.h"
//#include "IOPrep.h"
#include "VesselNetworkExporter.h"

typedef std::vector<Vessel> VesselVector;
typedef unsigned long ulong;

// initialize vascular network
Network::Network(int numLevels, int numDimensions)
	: numberOfLevels(numLevels)
	, dims(numDimensions)
{
	// set angle scaling constants for num_dimensions
	if (dims == 3)							
	{
		initialBifurcationAngle = 37.47;
		bifurcationReductionRatio = 0.88;
	}
	else if (dims == 2)
	{
		initialBifurcationAngle = 70;
		bifurcationReductionRatio = 0.6;
	}
	else
	{
		std::cout << "Must be 2 or 3 Dimensions" << std::endl;
		throw;
	}
	maxVesselRadius = minVesselRadius / pow(radiusReductionRatio, numberOfLevels - 1);
	maxVesselLength = minVesselLength / pow(lengthReductionRatio, numberOfLevels - 1);
}

// generate artery root sequentially
void Network::generate_artery_root(VesselVector& vessels, double x_projection, double y_projection)
{
	double x_extent = maxVesselLength;
	for(int i = 1; i < numberOfLevels; i++) {
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double length_factor = pow(lengthReductionRatio, i);
		x_extent += maxVesselLength * length_factor * cos(M_PI * bifAngle / 180);
	}
	double x_1 = - x_extent + x_projection;
	double y_1 = 0 + y_projection;

	vessels[0] = Vessel(0, Coordinate(- x_extent,0,0), Coordinate(- x_extent + x_projection, 0, 0), maxVesselLength, maxVesselRadius);
	vessels[0].add_nodes(0, 1);
}

// generate artery body sequentially
void Network::generate_artery_body(VesselVector& vessels, int level, double x_projection, double y_projection, double length, double radius)
{
	for (int i = pow(2, level - 1); i < pow(2, level); i++){
		int child1_idx = i * 2 - 1;
		int child2_idx = i * 2;
		Coordinate parent_end = vessels[i-1].get_endingPoint();
		Coordinate child1_end;
		Coordinate child2_end;

		if ((level % 2) || (getNumDimensions() == 2)) {
			child1_end = parent_end.add(x_projection, y_projection, 0);
			child2_end = parent_end.add(x_projection, -y_projection, 0);
		} else {
			child1_end = parent_end.add(x_projection, 0, y_projection);
			child2_end = parent_end.add(x_projection, 0, -y_projection);
		}
		vessels[child1_idx] = Vessel(child1_idx, parent_end, child1_end, length, radius);
		vessels[child1_idx].add_nodes(i, (i*2));
		vessels[child2_idx] = Vessel(child2_idx, parent_end, child2_end, length, radius);
		vessels[child2_idx].add_nodes(i, (i*2) + 1);
	}
}

// generate arterial tree sequentially
void Network::generate_artery(VesselVector& vessels) 
{
	generate_artery_root(vessels, maxVesselLength, 0);
	for (int i = 1; i < numberOfLevels; i++) {
		double length = maxVesselLength * pow(lengthReductionRatio, i);
		double radius = maxVesselRadius * pow(radiusReductionRatio, i);
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double x_projection = length * cos(M_PI * bifAngle / 180);
		double y_projection = length * sin(M_PI * bifAngle / 180);
		generate_artery_body(vessels, i, x_projection, y_projection, length, radius);
	}
}

// generate vein root sequentially
void Network::generate_vein_root(VesselVector& vessels, double x_projection, double y_projection)
{
	double x_extent = maxVesselLength;
	for(int i = 1; i < numberOfLevels; i++) {
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double length_factor = pow(lengthReductionRatio, i);
		x_extent += maxVesselLength * length_factor * cos(M_PI * bifAngle / 180);
	}
	int sink_vessel_idx = numNodes - 1;
	vessels[numVessels - 1] = Vessel(numVessels-1, Coordinate(x_extent,0,0), Coordinate(x_extent-x_projection,0,0), maxVesselLength, maxVesselRadius);
	//vessels[numVessels - 1] = Vessel(numVessels-1, Coordinate(x_extent-x_projection,0,0), Coordinate(x_extent,0,0), maxVesselLength, maxVesselRadius);
	vessels[numVessels - 1].add_nodes((sink_vessel_idx-1), sink_vessel_idx);
}

// generate vein body sequentially
void Network::generate_vein_body(VesselVector& vessels, int level, double x_projection, double y_projection, double length, double radius)
{
	for (int i = pow(2, level-1); i < pow(2, level); i++){
		int node_idx = numNodes - i - 1;
		int index = numVessels - i;
		int vessel1_idx = numVessels - i * 2 - 1;
		int vessel2_idx = numVessels - i * 2;
		Coordinate parent_end = vessels[index].get_endingPoint();
		Coordinate child1_end;
		Coordinate child2_end;

		if ((level % 2) || (getNumDimensions() == 2)) {
			child1_end = parent_end.add(-x_projection, y_projection, 0);
			child2_end = parent_end.add(-x_projection, -y_projection, 0);
		} else {
			child1_end = parent_end.add(-x_projection, 0, y_projection);
			child2_end = parent_end.add(-x_projection, 0, -y_projection);
		}
		vessels[vessel1_idx] = Vessel(vessel1_idx, parent_end, child1_end, length, radius);
		//vessels[vessel1_idx] = Vessel(vessel1_idx, child1_end, parent_end, length, radius);
		vessels[vessel1_idx].add_nodes((numNodes - i*2 - 2), node_idx);
		vessels[vessel2_idx] = Vessel(vessel2_idx, parent_end, child2_end, length, radius);
		//vessels[vessel2_idx] = Vessel(vessel2_idx, child2_end, parent_end, length, radius);
		vessels[vessel2_idx].add_nodes((numNodes - i*2 - 1), node_idx);
	}
}

// generate veinous tree sequentially
void Network::generate_vein(VesselVector& vessels) 
{
	generate_vein_root(vessels, maxVesselLength, 0);
	for (int i = 1; i < numberOfLevels; i++) {
		double length = maxVesselLength * pow(lengthReductionRatio, i);
		double radius = maxVesselRadius * pow(radiusReductionRatio, i);
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double x_projection = length * cos(M_PI * bifAngle / 180);
		double y_projection = length * sin(M_PI * bifAngle / 180);
		generate_vein_body(vessels, i, x_projection, y_projection, length, radius);
	}
}

/*
 * BEGIN PARALLEL FUNCTIONS
 *	- These functions generate the vessel network across n compute nodes. The root networks are created on 
 * 		every compute node, while each compute node creates its unique sub-networks
 * 	- Note: the vessel and junction node indexing is different than sequential creation due to
 * 		the nature of parallel generation
 */

// function to create the parallel root network - creates log2(mpiSize)+1 levels on each compute node
void Network::generate_root_parallel(VesselVector& arteries, VesselVector& veins, int rootNumLevels) {

	// calculate max x displacement
	double x_extent = maxVesselLength;
	for(int i = 1; i < numberOfLevels; i++) {
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double length_factor = pow(lengthReductionRatio, i);
		x_extent += maxVesselLength * length_factor * cos(M_PI * bifAngle / 180);
	}

	// generate initial vessels
	arteries.push_back(Vessel(Coordinate(- x_extent,0,0), Coordinate(- x_extent + maxVesselLength, 0, 0), maxVesselRadius, 0, 1));
	veins.push_back(Vessel(Coordinate(x_extent - maxVesselLength,0,0), Coordinate(x_extent,0,0), maxVesselRadius, numNodes - 2, numNodes - 1));
	// arteries.push_back(Vessel(Coordinate(- x_extent,0,0), Coordinate(- x_extent + maxVesselLength, 0, 0), maxVesselLength, maxVesselRadius));
	// arteries[0].add_nodes(0, 1);
	// veins.push_back(Vessel(Coordinate(x_extent - maxVesselLength,0,0), Coordinate(x_extent,0,0), maxVesselLength, maxVesselRadius));
	// veins[0].add_nodes(numNodes - 2, numNodes - 1);

	// calculate constants
	int subNumLevels = numberOfLevels - rootNumLevels;
	int rootNumVessels = pow(2, rootNumLevels-1) - 1;
	int subNumVessels = pow(2, subNumLevels+2) - 2;

	// loop over generation level to determine parameters
	for (int i = 1; i < rootNumLevels; i++) {
		double length = maxVesselLength * pow(lengthReductionRatio, i);
		double radius = maxVesselRadius * pow(radiusReductionRatio, i);
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
		double x_projection = length * cos(M_PI * bifAngle / 180);
		double y_projection = length * sin(M_PI * bifAngle / 180);

		// loops over parent vessels to create 2 daughters for each
		for (int j = 0; j < pow(2, i); j+=2) {
			long long index = pow(2, i) + j - 1;
			long long parent_index = index / 2;

			// calculate daughter artery parameters
			long long artery_startNode = arteries[parent_index].getNode2();
			//long long artery_startNode = arteries[parent_index].get_nodes()[1];
			long long artery1_endNode = 2 * artery_startNode;
			long long artery2_endNode = artery1_endNode + 1;
			Coordinate artery_parent_end = arteries[parent_index].get_endingPoint();
			Coordinate artery1_end;
			Coordinate artery2_end;

			// calculate daughter vein parameters
			long long vein_endNode = veins[parent_index].getNode1();
			//long long vein_endNode = veins[parent_index].get_nodes()[0];
			long long vein1_start_node = numNodes - pow(2, i + 1) + j;
			long long vein2_start_node = vein1_start_node + 1;
			Coordinate vein_parent_end = veins[parent_index].get_startingPoint();
			Coordinate vein1_end;
			Coordinate vein2_end;

			// check if 2D or 3D and set coords accordingly
			if ((i % 2) || (getNumDimensions() == 2)) {
				artery1_end = artery_parent_end.add(x_projection, y_projection, 0);
				artery2_end = artery_parent_end.add(x_projection, -y_projection, 0);
				vein1_end = vein_parent_end.add(-x_projection, y_projection, 0);
				vein2_end = vein_parent_end.add(-x_projection, -y_projection, 0);
			} else {
				artery1_end = artery_parent_end.add(x_projection, 0, y_projection);
				artery2_end = artery_parent_end.add(x_projection, 0, -y_projection);
				vein1_end = vein_parent_end.add(-x_projection, 0, y_projection);
				vein2_end = vein_parent_end.add(-x_projection, 0, -y_projection);
			}

			// calculate vessel indices
			long long arteryParentIndex, artery1_idx, artery2_idx, vein1_idx, vein2_idx, veinChildIndex;
			if (i < rootNumLevels-1) {
				arteryParentIndex = parent_index;
				artery1_idx = index;
				artery2_idx = index + 1;
				vein1_idx = 2 * rootNumVessels - pow(2, i + 1) + j + 1;
				veinChildIndex = 2 * rootNumVessels + 1 - ceil(static_cast<double>(2 * rootNumVessels - vein1_idx) / 2.0f);
                vein2_idx = vein1_idx + 1;
			} else {
				arteryParentIndex = parent_index;
				artery1_idx = 2 * rootNumVessels + (subNumVessels * (j));
                artery2_idx = artery1_idx + subNumVessels;
				vein1_idx = artery2_idx - 1;
                vein2_idx = vein1_idx + subNumVessels;
                veinChildIndex = 2 * rootNumVessels + 1 - ceil(static_cast<double>((pow(2, i + 1) - j - 1)) / 2.0f);
			}
			
			// set daughter arteries and add nodes
			arteries.push_back(Vessel(artery_parent_end, artery1_end, radius, artery_startNode, artery1_endNode));
			arteries.push_back(Vessel(artery_parent_end, artery2_end, radius, artery_startNode, artery2_endNode));
			// arteries.push_back(Vessel(artery_parent_end, artery1_end, length, radius));
			// arteries[index].add_nodes(artery_startNode, artery1_endNode);
			// arteries.push_back(Vessel(artery_parent_end, artery2_end, length, radius));
			// arteries[index+1].add_nodes(artery_startNode, artery2_endNode);

			// set daughter veins and add nodes
			veins.push_back(Vessel(vein1_end, vein_parent_end, radius, vein1_start_node, vein_endNode));
			veins.push_back(Vessel(vein2_end, vein_parent_end, radius, vein2_start_node, vein_endNode));
			// veins.push_back(Vessel(vein1_end, vein_parent_end, length, radius));
			// veins[index].add_nodes(vein1_start_node, vein_endNode);
			// veins.push_back(Vessel(vein2_end, vein_parent_end, length, radius));
			// veins[index+1].add_nodes(vein2_start_node, vein_endNode);

			//!! add connected vessels
			arteries[index].add_connectedVessel(arteryParentIndex);
            arteries[index].add_connectedVessel(artery2_idx);
            arteries[index + 1].add_connectedVessel(arteryParentIndex);
            arteries[index + 1].add_connectedVessel(artery1_idx);
            arteries[parent_index].add_connectedVessel(artery1_idx);
            arteries[parent_index].add_connectedVessel(artery2_idx);

            veins[index].add_connectedVessel(veinChildIndex);
            veins[index + 1].add_connectedVessel(veinChildIndex);
            veins[index].add_connectedVessel(vein2_idx);
            veins[index + 1].add_connectedVessel(vein1_idx);
            veins[parent_index].add_connectedVessel(vein1_idx);
            veins[parent_index].add_connectedVessel(vein2_idx);
		}
	}
}

void PruneNetworks(VesselVector& arteries, VesselVector& veins, int rootNumLevels)
{
    int mpiSize, mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    Vessel temp;
    ulong arteryIndex = pow(2, rootNumLevels - 1) - 1 + mpiRank;

    temp = arteries[arteryIndex];
    arteries.clear();
    arteries.push_back(temp);

    temp = veins[mpiRank];
	//temp = veins[arteryIndex];
	//temp.print_vessel();
    veins.clear();
    veins.push_back(temp);
}

// Function to create the vascular sub-networks on each node
void Network::generate_sub_parallel(VesselVector& arteries, VesselVector& veins, int rootNumLevels, int subNumLevels) {
	int mpiRank, mpiSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	PruneNetworks(arteries, veins, rootNumLevels);

	// define constants
	const long long rootNodes = pow(2, rootNumLevels);
	const long long newNodes = pow(2, subNumLevels) * 3 - 4;
	const long long subNumNodes = 3 * pow(2, subNumLevels+1) / 2;
	const long long subVesselsPerNode = pow(2, subNumLevels+1) - 1;

	long long vesselOffset = (pow(2, rootNumLevels) - 2) + mpiRank * (pow(2, subNumLevels+2)-2);
	long long nodeOffset = rootNodes + mpiRank * newNodes;

	// resize vectors
	arteries.resize(subVesselsPerNode);
	veins.resize(subVesselsPerNode);

	// loop over each level to calculate parameters
	for (int i = 1; i < subNumLevels+1; i++) {
		double length = minVesselLength / pow(lengthReductionRatio, subNumLevels-i);
		double radius = minVesselRadius / pow(radiusReductionRatio, subNumLevels-i);
		double bifAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, (rootNumLevels-1)+i);
		double x_projection = length * cos(M_PI * bifAngle / 180);
		double y_projection = length * sin(M_PI * bifAngle / 180);

		// loop over each branch to create 2 daughter vessels
		long long numVesselsInLevel = pow(2, i);
		for (int j = 0; j < numVesselsInLevel; j += 2) {
			long long index = pow(2, i) + j - 1;
			long long parent_index = index / 2;
			
			// calculate artery daughter information
			long long artery_startNode = arteries[parent_index].getNode2();
			//long long artery_startNode = arteries[parent_index].get_nodes()[1];
			long long artery1_endNode = nodeOffset + pow(2, i) + j - 2;
			long long artery2_endNode = artery1_endNode + 1;
			Coordinate artery_parent_end = arteries[parent_index].get_endingPoint();
			Coordinate artery1_end;
			Coordinate artery2_end;

			// calculate vein daughter information
			long long vein_endNode = veins[parent_index].getNode1();
			//long long vein_endNode = veins[parent_index].get_nodes()[0];
			long long vein1_startNode = subNumNodes - pow(2, i+1) + j + nodeOffset - 2;
			long long vein2_startNode = vein1_startNode + 1;
			Coordinate vein_parent_end = veins[parent_index].get_startingPoint();
			Coordinate vein1_end;
			Coordinate vein2_end;

			// determine if 2D or 3D and set vessel endpoints
			// 	- if 3D, alternate y & z axes
			if ((i % 2) || (getNumDimensions() == 2)) {
				artery1_end = artery_parent_end.add(x_projection, y_projection, 0);
				artery2_end = artery_parent_end.add(x_projection, -y_projection, 0);
				vein1_end = vein_parent_end.add(-x_projection, y_projection, 0);
				vein2_end = vein_parent_end.add(-x_projection, -y_projection, 0);
			} else {
				artery1_end = artery_parent_end.add(x_projection, 0, y_projection);
				artery2_end = artery_parent_end.add(x_projection, 0, -y_projection);
				vein1_end = vein_parent_end.add(-x_projection, 0, y_projection);
				vein2_end = vein_parent_end.add(-x_projection, 0, -y_projection);
			}

			// calculate vessel indices
			long long arteryParentIndex, artery1_idx, artery2_idx, vein1_idx, vein2_idx, veinChildIndex;
            arteryParentIndex = parent_index + vesselOffset;
            artery1_idx = index + vesselOffset;
            artery2_idx = artery1_idx + 1;
            vein1_idx = 2 * subVesselsPerNode - pow(2, i + 1) + j + 1;
            veinChildIndex = 2 * subVesselsPerNode + 1 - ceil(static_cast<double>(2 * subVesselsPerNode - vein1_idx) / 2.0f) + vesselOffset;
            vein1_idx += +vesselOffset;
            vein2_idx = vein1_idx + 1;

			// create the daughter arteries
			arteries[index] = Vessel(artery_parent_end, artery1_end, radius, artery_startNode, artery1_endNode);
			arteries[index+1] = Vessel(artery_parent_end, artery2_end, radius, artery_startNode, artery2_endNode);
			// arteries[index] = Vessel(artery_parent_end, artery1_end, length, radius);
			// arteries[index].add_nodes(artery_startNode, artery1_endNode);
			// arteries[index+1] = Vessel(artery_parent_end, artery2_end, length, radius);
			// arteries[index+1].add_nodes(artery_startNode, artery2_endNode);
			
			// create daughter veins
			veins[index] = Vessel(vein1_end, vein_parent_end, radius, vein1_startNode, vein_endNode);
			veins[index+1] = Vessel(vein2_end, vein_parent_end, radius, vein2_startNode, vein_endNode);
			// veins[index] = Vessel(vein1_end, vein_parent_end, length, radius);
			// veins[index].add_nodes(vein1_startNode, vein_endNode);
			// veins[index+1] = Vessel(vein2_end, vein_parent_end, length, radius);
			// veins[index+1].add_nodes(vein2_startNode, vein_endNode);

			//!! add connected vessels
			arteries[index].add_connectedVessel(arteryParentIndex);
            arteries[index].add_connectedVessel(artery2_idx);
            arteries[index + 1].add_connectedVessel(arteryParentIndex);
            arteries[index + 1].add_connectedVessel(artery1_idx);
            arteries[parent_index].add_connectedVessel(artery1_idx);
            arteries[parent_index].add_connectedVessel(artery2_idx);

            veins[index].add_connectedVessel(veinChildIndex);
            veins[index + 1].add_connectedVessel(veinChildIndex);
            veins[index].add_connectedVessel(vein2_idx);
            veins[index + 1].add_connectedVessel(vein1_idx);
            veins[parent_index].add_connectedVessel(vein1_idx);
            veins[parent_index].add_connectedVessel(vein2_idx);

			if (i == subNumLevels) {
                arteries[index].add_connectedVessel(vein1_idx);
                arteries[index + 1].add_connectedVessel(vein2_idx);
                veins[index].add_connectedVessel(artery1_idx);
                veins[index + 1].add_connectedVessel(artery2_idx);
            }
		}
	}
}

// generate the vascular network in parallel
void Network::generate_parallel(std::string filename) {
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	// calculate network constants
	int rootNumLevels = log2(mpiSize) + 1;
	int subNumLevels = numberOfLevels - rootNumLevels;
	int numVesselsPerArray = pow(2, rootNumLevels) - 1;
	int rootNumVessels = pow(2, rootNumLevels-1) - 1;

	if (!mpiRank) { 
		std::cout << "Root Num Levels = " << rootNumLevels << std::endl;
		std::cout << "Sub Num Levels = " << subNumLevels << std::endl;
	}

	// reserve memory for the vessel vectors
	VesselVector arteries, veins;
	arteries.reserve(numVesselsPerArray);
	veins.reserve(numVesselsPerArray);

	// create the root network
	double t1 = MPI_Wtime();
	generate_root_parallel(arteries, veins, rootNumLevels);
	if (!mpiRank) {
		std::cout << "-- Root Network Generation Time: " << (MPI_Wtime()-t1) << " seconds" << std::endl;
	}

	// sort veins by end node
	std::stable_sort(veins.begin(), veins.end(), CompareEndNodes2);
	
	// Export root vessels
	double t2 = MPI_Wtime();
	std::vector<double> bbox = getBoundingBox();
	hid_t fileId = RootVesselExporter(filename, arteries, veins, numberOfLevels, rootNumLevels, bbox);
	if (!mpiRank) {
		std::cout << "-- Root Network Export Time: " << (MPI_Wtime()-t2) << " seconds" << std::endl;
	}

	// create the sub network
	double t3 = MPI_Wtime();
	generate_sub_parallel(arteries, veins, rootNumLevels, subNumLevels);
	if (!mpiRank) {
		std::cout << "-- Sub Network Generation Time: " << (MPI_Wtime()-t3) << " seconds" << std::endl;
	}

	// sort sub veins by end node
	double sort2_time = MPI_Wtime();
	std::stable_sort(veins.begin(), veins.end(), CompareEndNodes2);
	if (!mpiRank) {
		std::cout << "-- Sub Network Sort Time: " << (MPI_Wtime()-sort2_time) << " seconds" << std::endl;
	}

	// export sub vessels
	double t4 = MPI_Wtime();
	SubVesselExporter(fileId, arteries, veins, numberOfLevels, subNumLevels);
	if (!mpiRank) {
		std::cout << "-- Sub Network Export Time: " << (MPI_Wtime()-t4) << " seconds" << std::endl;
	}

	// clear arteries, veins
	arteries.clear();
	veins.clear();
}

// Function call that generates the vascular network 
//	- automatically determines if application is running sequentially or 
// 		in parallel and builds accordingly
//? change function to void?
//VesselVector Network::generate(std::string filename)
void Network::generate(std::string filename)
{
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	// determine if generating sequentially or in parallel
	// VesselVector vessels;
	if (mpiSize == 1) {
		// reserve memory for vessel vector
		VesselVector vessels;
		vessels.resize(numVessels, Vessel());

		// generate arterial and venous trees
		double t1 = MPI_Wtime();
		generate_artery(vessels); 
	 	generate_vein(vessels);
		std::cout << "Sequential Network Generation Time: " << (MPI_Wtime()-t1) << " seconds" << std::endl;

		// export sequential
		double t2 = MPI_Wtime();
		std::vector<double> bbox = getBoundingBox();
		VesselNetworkExporter_HDF5(filename, numberOfLevels, vessels, bbox);
		std::cout << "Sequential Export Time: " << (MPI_Wtime()-t2) << " seconds" << std::endl;
	}
	else {
		// generate network in parallel
		generate_parallel(filename);
	}
	//return vessels;
}

// function to print the entire vector of vessels
void Network::print_vessels(VesselVector& vessels)
{
	for(int i = 0; i < getNumVessels(); i++){
		vessels[i].print_vessel();
	}
}

// function to calculate and return the output flow rate of the system
double Network::getOutputFlowRate()
{
	return pow(2, numberOfLevels-1) * M_PI * pow(minVesselRadius, 2);
}

std::vector<double> Network::getBoundingBox() {
	double xExtent = maxVesselLength;
	double yExtent = 0, zExtent = 0;
	double bifurcationAngle;
	double vesselLength = maxVesselLength;

	// determine if 2 or 3 dimensions and calculate extents accordingly
	if (getNumDimensions() == 2) {
		//zExtent = maxVesselRadius;
		for (int i = 1; i < numberOfLevels; i++) {
			bifurcationAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
			double power = numberOfLevels - 1 - i;
			vesselLength = minVesselLength / pow(lengthReductionRatio, power);

			xExtent += (vesselLength * cos(M_PI * bifurcationAngle / 180));
			yExtent += (vesselLength * sin(M_PI * bifurcationAngle / 180));
		}
	}
	else 
	{
		for (int i = 1; i < numberOfLevels; i++) {
        	bifurcationAngle = initialBifurcationAngle * pow(bifurcationReductionRatio, i);
			double power = numberOfLevels - 1 - i;
			vesselLength = minVesselLength / pow(lengthReductionRatio, power);
        	xExtent += (vesselLength * cos(M_PI * bifurcationAngle / 180));
        	if ((i % 2))
        	{
            	yExtent += (vesselLength * sin(M_PI * bifurcationAngle / 180));
        	}
        	else
        	{
            	zExtent += (vesselLength * sin(M_PI * bifurcationAngle / 180));
        	}
		}
	}
	//std::vector<double> bbox{xExtent, yExtent, zExtent};
	boundingBox = {-xExtent, xExtent, -yExtent, yExtent, -zExtent, zExtent};
	return boundingBox;
}

// class deconstructor
Network::~Network() {}