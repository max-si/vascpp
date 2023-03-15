#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include <mpi.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>

#include "BloodFlowEdgeNode.h"
#include "FlowCalculations.h"
#include "NetworkDescription.h"
#include "Triplet.h"

using std::cout;
using std::endl;
using std::vector;

typedef vector<BloodFlowVessel> VesselVector;
//typedef vector<Triplet> TripletList;
typedef unsigned long ulong;
typedef Eigen::Triplet<double> T;
//! changed all "Triplet" to "T"

/**
 * Function that creates the matrix elements for the flow rate
 *
 * Takes a member function pointer to specify what data the matrix elements will consist of.
 * @param[in] network
 * @param[in,out] tripletList
 * @param[in] edgeNodeMap
 */
void FlowListCreator(const NetworkDescription& network,
    std::vector<T>& tripletList,
    std::map<long long, long long>& edgeNodeMap,
    double (BloodFlowVessel::*conductanceTypePtr)() const)
{
    tripletList.resize(2 * network.getLocalNumberOfVessels());

    long long rowOffset = network.getLocalRankStartRow();
    const VesselVector& vesselVector = network.getLocalVesselVector();
    for (long long i = 0; i < vesselVector.size(); ++i)
    {
        long long index = i * 2;
        double conductance = (vesselVector[i].*conductanceTypePtr)();
        long long node = vesselVector[i].getNode1Index();
        if (vesselVector[i].getNode1IsLocal())
        {
            tripletList[index] = T(i, node - rowOffset, conductance);
        }
        else
        {
            long long newIndex = edgeNodeMap[node];
            tripletList[index] = T(i, newIndex, conductance);
        }
        node = vesselVector[i].getNode2Index();
        if (vesselVector[i].getNode2IsLocal())
        {
            tripletList[index + 1] = T(i, node - rowOffset, -conductance);
        }
        else
        {
            long long newIndex = edgeNodeMap[node];
            tripletList[index + 1] = T(i, newIndex, -conductance);
        }
    }
}

void HealthyFlowListCreator(const NetworkDescription& network,
    std::vector<T>& tripletList,
    std::map<long long, long long>& edgeNodeMap)
{
    FlowListCreator(network, tripletList, edgeNodeMap,
        &BloodFlowVessel::calculateHealthyConductance);
}

void DamagedFlowListCreator(const NetworkDescription& network,
    std::vector<T>& tripletList,
    std::map<long long, long long>& edgeNodeMap)
{
    FlowListCreator(network, tripletList, edgeNodeMap,
        &BloodFlowVessel::calculateDamagedConductance);
}

/**
 *  Function handles communication of pressures located on different compute nodes
 *
 *  Performs a "Gather All" opertaion of all the pressures for edge nodes. Add pressures to end of pressure vector.
 *  Updates the edge node map to map node index to correct location of pressure in the vector.
 *
 * @param[in] network
 * @param[in,out] pressureSolutionVector
 * @param[in] edgeNodeMap
 */
void GetRequiredPressuresFromOtherNodes(const NetworkDescription& network,
    std::vector<double>& pressureSolutionVector,
    std::map<long long, long long>& edgeNodeMap)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //determine which edge nodes are related to partition and sort into two sets, those to send and those to recieve
    vector<int> numNodePressureToRecievePerNode(mpiSize);
    vector<int> numNodePressureToSendPerNode(mpiSize);
    std::map<int, std::vector<int>> nodePressuresToSendMap;
    std::vector<std::pair<long long, double>> indexPressurePairsToSend;

    const std::set<BloodFlowEdgeNode>& edgeNodes = network.getEdgeNodeSet();
    for (auto itr = edgeNodes.begin(); itr != edgeNodes.end(); ++itr) {
        //cout << itr->getHostPartition() << " " << network.getLocalPartitionIndex();
        if (itr->getHostPartition() == network.getLocalPartitionIndex()) {
            auto partitions = itr->getAssociatedPartitions();
            vector<double>::size_type localPressureIndex = itr->getNodeIndex() - network.getLocalRankStartRow();
            indexPressurePairsToSend.emplace_back(itr->getNodeIndex(), pressureSolutionVector[localPressureIndex]);

            int pairToSendIndex = indexPressurePairsToSend.size() - 1;
            for (auto itr2 = partitions.begin(); itr2 != partitions.end(); itr2++) {
                if (nodePressuresToSendMap.count(*itr2)) {
                    (*nodePressuresToSendMap.find(*itr2)).second.push_back(pairToSendIndex);
                    numNodePressureToSendPerNode[*itr2] += 1;
                }
                else {
                    vector<int> temp;
                    temp.push_back(pairToSendIndex);
                    nodePressuresToSendMap.insert(std::make_pair(*itr2, temp));
                    numNodePressureToSendPerNode[*itr2] += 1;
                }
            }
        }
        if (itr->getAssociatedPartitions().count(network.getLocalPartitionIndex())) {
            numNodePressureToRecievePerNode[itr->getHostPartition()] += 1;
        }
    }
    int localTotalToSend = std::accumulate(numNodePressureToSendPerNode.begin(),
        numNodePressureToSendPerNode.end(), 0);
    int localTotalToReceive =
        std::accumulate(numNodePressureToRecievePerNode.begin(),
            numNodePressureToRecievePerNode.end(), 0);
    //construct information for messages
    vector<long long> indexesToSend(localTotalToSend);
    vector<double> pressuresToSend(localTotalToSend);

    vector<long long> recievedIndexes(localTotalToReceive);
    vector<double> receivedPressures(localTotalToReceive);

    vector<int> sendOffsets;
    sendOffsets.push_back(0);
    vector<int> recvOffsets;
    recvOffsets.push_back(0);

    std::partial_sum(numNodePressureToSendPerNode.begin(),
        numNodePressureToSendPerNode.end() - 1,
        std::back_inserter(sendOffsets));
    std::partial_sum(numNodePressureToRecievePerNode.begin(),
        numNodePressureToRecievePerNode.end() - 1,
        std::back_inserter(recvOffsets));
    int index = 0;
    for (std::map<int, vector<int>>::iterator itr = nodePressuresToSendMap.begin(); itr != nodePressuresToSendMap.end(); ++itr) {
        const vector<int>& sendIndexes = itr->second;
        for (int i = 0; i < sendIndexes.size(); ++i) {
            indexesToSend[index] = indexPressurePairsToSend[sendIndexes[i]].first;
            pressuresToSend[index] = indexPressurePairsToSend[sendIndexes[i]].second;
            index++;
        }
    }
    //Communication 
    MPI_Alltoallv(indexesToSend.data(), numNodePressureToSendPerNode.data(),
        sendOffsets.data(), MPI_LONG_LONG_INT, recievedIndexes.data(),
        numNodePressureToRecievePerNode.data(), recvOffsets.data(),
        MPI_LONG_LONG_INT, MPI_COMM_WORLD);
    MPI_Alltoallv(pressuresToSend.data(), numNodePressureToSendPerNode.data(),
        sendOffsets.data(), MPI_DOUBLE, receivedPressures.data(),
        numNodePressureToRecievePerNode.data(), recvOffsets.data(), MPI_DOUBLE,
        MPI_COMM_WORLD);
    //append node num to pressure list

    long long numLocalPressureRows = network.getLocalRankEndRow() - network.getLocalRankStartRow() + 1;
    for (int i = 0; i < localTotalToReceive; ++i)
    {
        edgeNodeMap.insert(std::make_pair(recievedIndexes[i], numLocalPressureRows + i));
        pressureSolutionVector.push_back(receivedPressures[i]);
    }
}


void HealthyFlowCalculations(NetworkDescription& network, std::vector<double>& pressureSolutionVector)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    vector<T> flowTripletList;
    std::map<long long, long long> edgeNodeMap;
    GetRequiredPressuresFromOtherNodes(network, pressureSolutionVector, edgeNodeMap);
    HealthyFlowListCreator(network, flowTripletList, edgeNodeMap);

    // for (auto& x : flowTripletList)
	// {
	// 	cout << x << '\n';
	// } cout << endl;

    // for(int i = 0; i < pressureSolutionVector.size(); i++) {
	// 	std::cout << pressureSolutionVector[i] << std::endl;
	// }

    // solve for flows of each vessel
    Eigen::SparseMatrix<double> PoiseuilleMatrix(network.getLocalNumberOfVessels(), pressureSolutionVector.size());
	PoiseuilleMatrix.setFromTriplets(flowTripletList.begin(), flowTripletList.end());
    flowTripletList.clear();

    Eigen::VectorXd pressureVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(pressureSolutionVector.data(), pressureSolutionVector.size());
    Eigen::VectorXd flowVector = PoiseuilleMatrix * pressureVector;
    PoiseuilleMatrix.resize(0,0);
    pressureVector.resize(0);

    std::vector<double> flowSolutionVector(flowVector.data(), flowVector.data() + flowVector.size());

    // for(int i = 0; i < flowSolutionVector.size(); i++) {
	// 	std::cout << flowSolutionVector[i] << std::endl;
	// }
    // std::cout << mpiRank << ": " << flowSolutionVector.size() << std::endl;

	
	
	// for (int i = 0; i < mpiSize; i++)
	// {
	// 	if (i == mpiRank)
	// 	{
	// 		cout << "Rank " << i << " Blood Flows:\n";
	// 		for (size_t j = 0; j < network.getLocalNumberOfVessels(); j++)
	// 			cout << flowSolutionVector[j] << "\n";
	// 			cout <<endl;
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
    
    AddHealthyFlowToVessels(network.getLocalVesselVector(), flowSolutionVector);
}


/**
 * \brief Function that adds flow rates to vessel objects
 *
 * Takes advantage of function pointers to allow a single function to perform both duties.
 * Put into wrapper functions to make coding easier.
 * In future could be a template function for better resuability.
 *
 * @param[in,out] vessels
 * @param[in] flowValues
 * @param[in] flowTypePtr
 */
void AddFlowValuesToVessels(std::vector<BloodFlowVessel>& vessels,
    std::vector<double>& flowValues,
    void (BloodFlowVessel::*flowTypePtr)(double));

void AddHealthyFlowToVessels(
    std::vector<BloodFlowVessel>& vessels, std::vector<double>& flowValues)
{
    AddFlowValuesToVessels(
        vessels, flowValues, &BloodFlowVessel::SetHealthyFlow);
}

void AddDamagedFlowToVessels(
    std::vector<BloodFlowVessel>& vessels, std::vector<double>& flowValues)
{
    AddFlowValuesToVessels(
        vessels, flowValues, &BloodFlowVessel::SetDamagedFlow);
}

void AddFlowValuesToVessels(std::vector<BloodFlowVessel>& vessels,
    std::vector<double>& flowValues,
    void (BloodFlowVessel::*flowTypePtr)(double))
{
    ulong index = 0;
    VesselVector::iterator itr;

    for (itr = vessels.begin(); itr != vessels.end(); itr++)
    {
        ((*itr).*flowTypePtr)(flowValues[index++]);
    }
}