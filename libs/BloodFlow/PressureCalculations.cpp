#include <iostream>
#include <utility>
#include <vector>

#include <mpi.h>

#include "BloodFlowVessel.h"
#include "NetworkDescription.h"
#include "PressureCalculations.h"
#include "Triplet.h"
#include "PressureListCreator.h"
#include "Solver.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef vector<BloodFlowVessel> VesselVector;
typedef vector<Triplet> TripletList;

void CalculatePressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector,
    void (*PressureListFunction)(const NetworkDescription&, TripletList&));

void CalculateDamagedPressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector)
{
    return CalculatePressures(network, pressureSolutionVector, DamagedPressureListCreator);
}

void CalculateHealthyPressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector)
{
    return CalculatePressures(network, pressureSolutionVector, HealthyPressureListCreator);
}

/**
 * Generic function that implementes matrix solving pipeline to find pressures in network junctions
 *
 * Means changes to this function apply to both the healthy and damaged calculations
 *
 * @param network
 * @param pressureSolutionVector
 * @param PressureListFunction -takes a pressure list function pointer to create pressure lists
 */
void CalculatePressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector,
    void (*PressureListFunction)(const NetworkDescription&, TripletList&))
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    TripletList tripletList;

    MPI_Barrier(MPI_COMM_WORLD);
    //if (!mpiRank) { cout << "PressureListCreator" << endl; }
    PressureListFunction(network, tripletList);

    //if (!mpiRank) { cout << "Permuting Matrix" << endl; }
    MPI_Barrier(MPI_COMM_WORLD);
    CreatePressureRhsVector(network, pressureSolutionVector);

    if (!mpiRank) { cout << endl << "Entering Trilinos Solver" << endl; }
    MPI_Barrier(MPI_COMM_WORLD);
    SolvePressures(network, tripletList, pressureSolutionVector);

    MPI_Barrier(MPI_COMM_WORLD);
}

// create the pressure boundary conditions vector
void CreatePressureRhsVector(const NetworkDescription& network, std::vector<double>& rhsSolutionVector)
{
    long long numElementsInVector = network.getLocalRankEndRow() - network.getLocalRankStartRow() + 1;
    rhsSolutionVector.assign(numElementsInVector, 0.0f);

    int index = 0;

    if (!network.GetPermutedBcNodeNumbers().size())
    {
        for (auto& x : network.getBoundaryConditionMap())
        {
            if (x.first >= network.getLocalRankStartRow() && x.first <= network.getLocalRankEndRow())
            {
                long long index = x.first - network.getLocalRankStartRow();
                rhsSolutionVector[index] = x.second.second;
            }
        }
    }
    else
    {
        for (auto& x : network.getBoundaryConditionMap())
        {
            long long nodeIndex = network.GetPermutedBcNodeNumbers()[index++];
            if (nodeIndex >= network.getLocalRankStartRow() && nodeIndex <= network.getLocalRankEndRow())
            {
                long long index = nodeIndex - network.getLocalRankStartRow();
                rhsSolutionVector[index] = x.second.second;
            }
        }
    }
}