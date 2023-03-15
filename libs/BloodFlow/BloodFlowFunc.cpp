

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
#include "FlowCalculations.h"
#include "BloodFlowFunc.h"
#include "BloodFlowExporter.h"

typedef std::vector<Triplet> TripletList;

void BloodFlowExe(std::string filename) {
	int mpiRank;
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	//import metadata
	NetworkDescription vesselNetwork;

	MPI_Barrier(MPI_COMM_WORLD);

	//if (!mpiRank) { std::cout << "Entering Importer" << std::endl; }
	hid_t fileId = BloodFlowNetworkImporter(filename, vesselNetwork);

	//std::cout << vesselNetwork << std::endl;

	// if (damageIndex > -1) {
	// 	DamageSpecificVesselAllSegment(vesselNetwork, damageIndex, 0.0);
	// }

	std::vector<double> pressureSolutionVector;
	CalculateHealthyPressures(vesselNetwork, pressureSolutionVector);
	MPI_Barrier(MPI_COMM_WORLD);
	// std::cout << mpiRank << ": " << pressureSolutionVector.size() << std::endl;

	if (!mpiRank) { std::cout << std::endl << "Calculating Healthy Flows" << std::endl; }
	HealthyFlowCalculations(vesselNetwork, pressureSolutionVector);

	MPI_Barrier(MPI_COMM_WORLD);

	// check if healthy flows are correct
	bool healthyFlowsCorrect = CheckHealthyFlowsAreCorrect(vesselNetwork);
	if (healthyFlowsCorrect) {
        std::cout << mpiRank << ": All healthy Flows calculated correctly" << std::endl;
    } else {
        std::cout << mpiRank << ": Healthy flows not inverse of power of two" << std::endl;
    }

	// calculate damaged pressures and flows
	
	// export flows to HDF5
	MPI_Barrier(MPI_COMM_WORLD);
	BloodFlowExporter(fileId, vesselNetwork);


	H5Fclose(fileId);
}


bool CheckHealthyFlowsAreCorrect(const NetworkDescription& network) {
	const std::vector<BloodFlowVessel>& vessels = network.getLocalVesselVector();

	//find flow boundary condition
	double bcFlowRate;
	for (auto& x : network.getBoundaryConditionMap()) {
		if (x.second.first == NetworkDescription::BoundrayConditionType::Flow)
			bcFlowRate = x.second.second;
	}

	bool isCorrect = true;
	for (long long i = 0; i < vessels.size(); ++i) {
		double power = bcFlowRate / vessels[i].GetHealthyFlow();
		if (fabs(power - int(std::round(power))) / int(std::round(power)) >	0.005) {
			std::cout << i << ": " << std::scientific
					<< vessels[i].GetHealthyFlow() << " " << power - int(power)
					<< " " << power << " " << int(std::round(power)) << " "
					<< std::endl;
			return false;
		}
		//cout <<std::scientific <<mantissa << endl;
	}

	return isCorrect;
}