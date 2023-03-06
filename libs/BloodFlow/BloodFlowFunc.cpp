

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "Triplet.h"
#include "NetworkDescription.h"
#include "NetworkImporter.h"
#include "PressureCalculations.h"

typedef std::vector<Triplet> TripletList;


void BloodFlowExe(std::string filename) {
	int mpiRank;
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	//import metadata
	NetworkDescription vesselNetwork;

	MPI_Barrier(MPI_COMM_WORLD);

	if (!mpiRank) { std::cout << "Entering Importer" << std::endl; }
	hid_t fileId = BloodFlowNetworkImporter(filename, vesselNetwork);

	std::cout << vesselNetwork << std::endl;

	// if (damageIndex > -1) {
	// 	DamageSpecificVesselAllSegment(vesselNetwork, damageIndex, 0.0);
	// }

	std::vector<double> pressureSolutionVector;
	CalculateHealthyPressures(vesselNetwork, pressureSolutionVector);
	MPI_Barrier(MPI_COMM_WORLD);

	//HealthyFlowCalculations(vesselNetwork, pressureSolutionVector);
	MPI_Barrier(MPI_COMM_WORLD);

	// check if healthy flows are correct

	// calculate damaged pressures and flows
	// MPI_Barrier(MPI_COMM_WORLD);

	// blood flow exporter



	H5Fclose(fileId);
}