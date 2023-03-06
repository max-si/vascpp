#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include <mpi.h>
//#include <omp.h>

#include "BloodFlowVessel.h"
#include "NetworkDescription.h"
#include "PressureListCreator.h"
#include "Triplet.h"
#include "core/Timsort.h"
#include "MergeTriplets.h"
// #include <vascular/core/ParallelSortTest.h>

using std::cout;
using std::endl;
using std::vector;

typedef vector<BloodFlowVessel> VesselVector;
typedef vector<Triplet> TripletArray;
typedef vector<Triplet> TripletList;
//typedef unsigned long ulong;

void CreatePressureTripletList(const NetworkDescription& network, TripletArray& localTripletArray, TripletArray& nonLocalTripletArray,
    double (BloodFlowVessel::*conductanceTypePtr)() const);
void PressureListCreator(const NetworkDescription& network, std::vector<Triplet>& tripletList,
    double (BloodFlowVessel::*conductanceTypePtr)() const);

void HealthyPressureListCreator(const NetworkDescription& network, std::vector<Triplet>& tripletList)
{
    PressureListCreator(network, tripletList, &BloodFlowVessel::calculateHealthyConductance);
}

void DamagedPressureListCreator(const NetworkDescription& network, std::vector<Triplet>& tripletList)
{
    PressureListCreator(network, tripletList, &BloodFlowVessel::calculateDamagedConductance);
}

/**
 * Function to communicate the non-local triplets in the network for creation of matrix elements
 *
 * Uses two MPI_AlltoAllv Calls to perform communication.
 * Could be simplified but this allows the most generic implementation
 * Automatically returns if there is only one partition.
 * @param network
 * @param nonLocalTriplets
 */
void NonLocalTripletCommunication(
    const NetworkDescription& network, std::vector<Triplet>& nonLocalTriplets)
{
    if (network.getNumberOfPartitions() == 1)
        return;

    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    vector<int> numTripletsToSendToEachNode;
    vector<int> numTripletsToRecieveFromEachNode;
    numTripletsToSendToEachNode.assign(mpiSize, 0);
    numTripletsToRecieveFromEachNode.assign(mpiSize, 0);

    std::map<int, std::vector<Triplet>> nodeTripletMap;

    // how many triplets do i need to send to other nodes
    long localTotalNumberToSend = nonLocalTriplets.size();
    for (std::vector<Triplet>::iterator itr = nonLocalTriplets.begin(); itr != nonLocalTriplets.end(); ++itr)
    {
        long long nodeIndex = itr->first;
        int hostPartitionIndex = network.getEdgeNodeHostPartition(nodeIndex);
        numTripletsToSendToEachNode[hostPartitionIndex] += 1;
    }

    vector<int> sendMatrix;
    sendMatrix.assign(mpiSize * mpiSize, 0);

    MPI_Allgather(numTripletsToSendToEachNode.data(), mpiSize, MPI_INT, sendMatrix.data(), mpiSize, MPI_INT, MPI_COMM_WORLD);
    long localTotalToRecieve = 0;

    for (int i = 0; i < mpiSize; ++i)
    {
        localTotalToRecieve += sendMatrix[i * mpiSize + mpiRank];
        numTripletsToRecieveFromEachNode[i] = sendMatrix[i * mpiSize + mpiRank];
    }

    vector<long long> indexesToSend(localTotalNumberToSend * 2);
    vector<long long> recievedIndexes(localTotalToRecieve * 2);
    vector<double> valuesToSend(localTotalNumberToSend);
    vector<double> valuesToRecieve(localTotalToRecieve);
    vector<int> valueSendOffsets;
    valueSendOffsets.push_back(0);
    vector<int> valueReceiveOffsets;
    valueReceiveOffsets.push_back(0);

    std::partial_sum(numTripletsToSendToEachNode.begin(),
        numTripletsToSendToEachNode.end() - 1,
        std::back_inserter(valueSendOffsets));
    std::partial_sum(numTripletsToRecieveFromEachNode.begin(),
        numTripletsToRecieveFromEachNode.end() - 1,
        std::back_inserter(valueReceiveOffsets));

    const vector<int>& valueSendOffsetRef = valueSendOffsets;
    const vector<int>& valueReceiveOffsetRef = valueReceiveOffsets;

    vector<int> indexSendOffsets(valueSendOffsetRef);
    vector<int> indexReceiveOffsets(valueReceiveOffsetRef);
    vector<int> numIndexesToSendPerProcess(numTripletsToSendToEachNode);
    vector<int> numIndexesToRecievePerProcess(numTripletsToRecieveFromEachNode);
    for (int i = 0; i < mpiSize; ++i)
    {
        indexReceiveOffsets[i] *= 2;
        indexSendOffsets[i] *= 2;
        numIndexesToSendPerProcess[i] *= 2;
        numIndexesToRecievePerProcess[i] *= 2;
    }

    int nodeIndex = 0, valueIndex = 0;
   
    for (vector<Triplet>::const_iterator itr = nonLocalTriplets.begin();
         itr != nonLocalTriplets.end(); ++itr)
    {
        indexesToSend[nodeIndex++] = itr->first;
        indexesToSend[nodeIndex++] = itr->second;
        valuesToSend[valueIndex++] = itr->third;
    }

    //cout << "Communicating Indexes" << endl;
    MPI_Alltoallv(indexesToSend.data(), numIndexesToSendPerProcess.data(),
        indexSendOffsets.data(), MPI_LONG_LONG_INT, recievedIndexes.data(),
        numIndexesToRecievePerProcess.data(), indexReceiveOffsets.data(),
        MPI_LONG_LONG_INT, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout << "Communicating Values" << endl;

    MPI_Alltoallv(valuesToSend.data(), numTripletsToSendToEachNode.data(),
        valueSendOffsets.data(), MPI_DOUBLE, valuesToRecieve.data(),
        numTripletsToRecieveFromEachNode.data(), valueReceiveOffsets.data(),
        MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    nonLocalTriplets.clear();
    for (int i = 0; i < localTotalToRecieve; ++i)
    {
        int nodeIndex2 = i * 2;
        long long node1 = recievedIndexes[nodeIndex2];
        long long node2 = recievedIndexes[nodeIndex2 + 1];
        double value = valuesToRecieve[i];
        nonLocalTriplets.push_back(Triplet(node1, node2, value));
    }

    gfx::timsort(nonLocalTriplets.begin(), nonLocalTriplets.end(), CompareTriplet);
    MergeTripletsInPlace(nonLocalTriplets);
}

#include <set>
void TestTripletListCreator(const TripletArray& localTripletArray, const NetworkDescription& network)
{
    std::multiset<Triplet> testSet;

    const vector<BloodFlowVessel>& vessels = network.getLocalVesselVector();

    for (vector<BloodFlowVessel>::const_iterator itr = vessels.begin();
         itr != vessels.end(); ++itr)
    {
        Node node1 = itr->getNode1();
        Node node2 = itr->getNode2();
        double conductance = itr->calculateHealthyConductance();
        //node 1
        if (node1.GetIsBcNode())
        {
            const std::pair<NetworkDescription::BoundrayConditionType, double>
                bcInfo = network.getBoundaryConditionNodeInformation(node1.GetIndex());
            if (bcInfo.first == NetworkDescription::BoundrayConditionType::Pressure)
            {
                Triplet temp(node1.GetIndex(), node1.GetIndex(), 1);
                testSet.insert(temp);
                temp = Triplet(node1.GetIndex(), node2.GetIndex(), 0);
                testSet.insert(temp);
            }
            else
            {
                Triplet temp(node1.GetIndex(), node1.GetIndex(), -conductance);
                testSet.insert(temp);
                temp = Triplet(node1.GetIndex(), node2.GetIndex(), conductance);
                testSet.insert(temp);
            }
        }
        else
        {
            Triplet temp(node1.GetIndex(), node1.GetIndex(), -conductance);
            testSet.insert(temp);
            temp = Triplet(node1.GetIndex(), node2.GetIndex(), conductance);
            testSet.insert(temp);
        }
        //node 2
        if (node2.GetIsBcNode())
        {
            const std::pair<NetworkDescription::BoundrayConditionType, double>
                bcInfo = network.getBoundaryConditionNodeInformation(node2.GetIndex());
            if (bcInfo.first == NetworkDescription::BoundrayConditionType::Pressure)
            {
                Triplet temp(node2.GetIndex(), node2.GetIndex(), 1);
                testSet.insert(temp);
                temp = Triplet(node2.GetIndex(), node1.GetIndex(), 0);
                testSet.insert(temp);
            }
            else
            {
                Triplet temp(node2.GetIndex(), node2.GetIndex(), -conductance);
                testSet.insert(temp);
                temp = Triplet(node2.GetIndex(), node1.GetIndex(), conductance);
                testSet.insert(temp);
            }
        }
        else
        {
            Triplet temp(node2.GetIndex(), node2.GetIndex(), -conductance);
            testSet.insert(temp);
            temp = Triplet(node2.GetIndex(), node1.GetIndex(), conductance);
            testSet.insert(temp);
        }
    }

    // check if the vectors match
    if (std::equal(testSet.begin(), testSet.end(), localTripletArray.begin()))
    {
        cout << "triplet lists are equivalent" << endl;
    }

    else
    {
        cout << "Size: " << testSet.size() << " " << localTripletArray.size() << endl;
        cout << "triplet lists are not equivalent" << endl;

        long long errorCount = 0;
        vector<Triplet>::const_iterator itr = localTripletArray.begin();
        for (std::multiset<Triplet>::const_iterator setItr = testSet.begin(); setItr != testSet.end(); ++setItr, ++itr)
        {
            errorCount += (*setItr != *itr);
            //cout << *setItr <<"\t" <<*itr<<endl;
        }
        cout << errorCount << endl;
        throw;
    }
}

void TestMergedTripletListCreator(
    const TripletArray& localTripletArray, const NetworkDescription& network)
{
    std::set<Triplet> testSet;

    const vector<BloodFlowVessel>& vessels = network.getLocalVesselVector();

    for (vector<BloodFlowVessel>::const_iterator itr = vessels.begin(); itr != vessels.end(); ++itr)
    {
        Node node1 = itr->getNode1();
        Node node2 = itr->getNode2();
        double conductance = itr->calculateHealthyConductance();
        //node 1
        if (node1.GetIsBcNode())
        {
            const std::pair<NetworkDescription::BoundrayConditionType, double>
                bcInfo = network.getBoundaryConditionNodeInformation(node1.GetIndex());
            if (bcInfo.first == NetworkDescription::BoundrayConditionType::Pressure)
            {
                Triplet temp(node1.GetIndex(), node1.GetIndex(), 1);
                if (testSet.count(temp))
                {
                    auto tripletItr = testSet.find(temp);
                    temp += *tripletItr;
                    testSet.erase(tripletItr);
                    testSet.insert(temp);
                }
                else
                {
                    testSet.insert(temp);
                }

                temp = Triplet(node1.GetIndex(), node2.GetIndex(), 0);
                testSet.insert(temp);
            }
            else
            {
                Triplet temp(node1.GetIndex(), node1.GetIndex(), -conductance);
                if (testSet.count(temp))
                {
                    auto tripletItr = testSet.find(temp);
                    temp += *tripletItr;
                    testSet.erase(tripletItr);
                    testSet.insert(temp);
                }
                else
                {
                    testSet.insert(temp);
                }
                temp = Triplet(node1.GetIndex(), node2.GetIndex(), conductance);
                testSet.insert(temp);
            }
        }
        else
        {
            Triplet temp(node1.GetIndex(), node1.GetIndex(), -conductance);
            if (testSet.count(temp))
            {
                auto tripletItr = testSet.find(temp);
                temp += *tripletItr;
                testSet.erase(tripletItr);
                testSet.insert(temp);
            }
            else
            {
                testSet.insert(temp);
            }
            temp = Triplet(node1.GetIndex(), node2.GetIndex(), conductance);
            testSet.insert(temp);
        }
        //node 2
        if (node2.GetIsBcNode())
        {
            const std::pair<NetworkDescription::BoundrayConditionType, double>
                bcInfo = network.getBoundaryConditionNodeInformation(node2.GetIndex());
            if (bcInfo.first == NetworkDescription::BoundrayConditionType::Pressure)
            {
                Triplet temp(node2.GetIndex(), node2.GetIndex(), 1);
                if (testSet.count(temp))
                {
                    auto tripletItr = testSet.find(temp);
                    temp += *tripletItr;
                    testSet.erase(tripletItr);
                    testSet.insert(temp);
                }
                else
                {
                    testSet.insert(temp);
                }
                temp = Triplet(node2.GetIndex(), node1.GetIndex(), 0);
                testSet.insert(temp);
            }
            else
            {
                Triplet temp(node2.GetIndex(), node2.GetIndex(), -conductance);
                if (testSet.count(temp))
                {
                    auto tripletItr = testSet.find(temp);
                    temp += *tripletItr;
                    testSet.erase(tripletItr);
                    testSet.insert(temp);
                }
                else
                {
                    testSet.insert(temp);
                }
                temp = Triplet(node2.GetIndex(), node1.GetIndex(), conductance);
                testSet.insert(temp);
            }
        }
        else
        {
            Triplet temp(node2.GetIndex(), node2.GetIndex(), -conductance);
            if (testSet.count(temp))
            {
                auto tripletItr = testSet.find(temp);
                temp += *tripletItr;
                testSet.erase(tripletItr);
                testSet.insert(temp);
            }
            else
            {
                testSet.insert(temp);
            }
            temp = Triplet(node2.GetIndex(), node1.GetIndex(), conductance);
            testSet.insert(temp);
        }
    }

    // check if the vectors match

    if (std::equal(testSet.begin(), testSet.end(), localTripletArray.begin()))
    {
        cout << "triplet lists are equivalent" << endl;
    }

    else
    {
        cout << "Size: " << testSet.size() << " " << localTripletArray.size() << endl;
        cout << "triplet lists are not equivalent" << endl;

        long long errorCount = 0;
        vector<Triplet>::const_iterator itr = localTripletArray.begin();
        for (std::multiset<Triplet>::const_iterator setItr = testSet.begin(); setItr != testSet.end(); ++setItr, ++itr)
        {
            errorCount += (*setItr != *itr);
            //cout << *setItr <<"\t" <<*itr<<endl;
        }
        cout << errorCount << endl;
        throw;
    }
}

/**
 * Function that creates the pressure triplets for the matrix elements
 *
 * Function takes a member function pointer specify which data to use to create the triplet
 * @param network
 * @param tripletList
 * @param conductanceTypePtr
 */
void PressureListCreator(const NetworkDescription& network,
    std::vector<Triplet>& tripletList,
    double (BloodFlowVessel::*conductanceTypePtr)() const)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    tripletList.resize(4 * (network.getLocalNumberOfVessels() + 2));

    TripletArray& localTripletArray = tripletList;
    TripletArray nonLocalTripletArray;
    //localTripletArray.resize(network.getLocalNumberOfVessels()* 4+2);
    
    CreatePressureTripletList(network, localTripletArray, nonLocalTripletArray, conductanceTypePtr);
    /*for (auto& x : localTripletArray)
	{
		cout << x << '\n';
	}
	cout <<endl;*/
    //cout << "nonLocalTripletArray size " <<nonLocalTripletArray.size()<<endl;


    gfx::timsort(nonLocalTripletArray.begin(), nonLocalTripletArray.end(), CompareTriplet);
    MergeTripletsInPlace(nonLocalTripletArray);


// #pragma omp parallel shared(localTripletArray)
//     {
// #pragma omp master
//         {
//             pss::parallel_stable_sort(localTripletArray.begin(),
//                 localTripletArray.end(), CompareTriplet);
//         }
//     }

    /*for (auto& x : localTripletArray)
	{
		cout << x << '\n';
	}
	cout << endl;*/

    Triplet emptyTriplet = Triplet();
    TripletArray::iterator onePastEmptyTriplets = std::upper_bound(localTripletArray.begin(), localTripletArray.end(), emptyTriplet, CompareTriplet);
    localTripletArray.erase(localTripletArray.begin(), onePastEmptyTriplets);
    //cout << "after delete" <<endl;
    /*for (auto& x : localTripletArray)
	{
		cout << x << '\n';
	}*/

    //TestTripletListCreator(localTripletArray,network);

    MergeTriplets(localTripletArray);

    // communicate non-local nodes to where they belong
    NonLocalTripletCommunication(network, nonLocalTripletArray);

    for (vector<Triplet>::iterator itr = nonLocalTripletArray.begin();
         itr != nonLocalTripletArray.end(); ++itr)
    {
        if (itr->first == itr->second)    // A host partition will always have a diagonal element
        {
            vector<Triplet>::iterator tripletListItr = std::find(tripletList.begin(), tripletList.end(), *itr);
            tripletListItr->third += itr->third;
        }
        else
        {
            tripletList.push_back(*itr);
        }
    }

    gfx::timsort(tripletList.begin(), tripletList.end(), CompareTriplet);

    for (TripletArray::iterator itr = tripletList.begin(); itr != tripletList.end(); ++itr)
    {
        itr->first -= network.getLocalRankStartRow();
    }
}

void CreatePressureTripletList(const NetworkDescription& network,
    TripletArray& localTripletArray, TripletArray& nonLocalTripletArray,
    double (BloodFlowVessel::*conductanceTypePtr)() const)
{
    const VesselVector& vessels = network.getLocalVesselVector();

    for (long long i = 0; i < vessels.size(); ++i)
    {
        long long index = i * 4;
        double conductance = ((vessels[i]).*conductanceTypePtr)();
        const long long node1Index = vessels[i].getNode1Index();
        const long long node2Index = vessels[i].getNode2Index();
        //check node 1
        if (vessels[i].getNode1IsBcNode() &&
            network.getBoundaryConditionNodeInformation(node1Index).first == NetworkDescription::Pressure)
        {
            if (vessels[i].getNode1IsLocal())
            {
                localTripletArray[index] = Triplet(node1Index, node1Index, 1.0f);
                localTripletArray[index + 1] = Triplet(node1Index, node2Index, 0.0f);
            }
            else
            {
                nonLocalTripletArray.emplace_back(node1Index, node1Index, 1.0f);
                nonLocalTripletArray.emplace_back(node1Index, node2Index, 0.0f);
            }
        }
        else
        {
            if (vessels[i].getNode1IsLocal())
            {
                localTripletArray[index] = Triplet(node1Index, node1Index, -conductance);
                localTripletArray[index + 1] = Triplet(node1Index, node2Index, conductance);
            }
            else
            {
                nonLocalTripletArray.emplace_back(node1Index, node1Index, -conductance);
                nonLocalTripletArray.emplace_back(node1Index, node2Index, conductance);
            }
        }
        //check node 2
        if (vessels[i].getNode2IsBcNode() &&
            network.getBoundaryConditionNodeInformation(node2Index).first == NetworkDescription::Pressure)
        {
            if (vessels[i].getNode2IsLocal())
            {
                localTripletArray[index + 2] = Triplet(node2Index, node2Index, 1.0f);
                localTripletArray[index + 3] = Triplet(node2Index, node1Index, 0.0f);
            }
            else
            {
                nonLocalTripletArray.emplace_back(node2Index, node2Index, 1.0f);
                nonLocalTripletArray.emplace_back(node2Index, node1Index, 0.0f);
            }
        }
        else
        {
            if (vessels[i].getNode2IsLocal())
            {
                localTripletArray[index + 2] = Triplet(node2Index, node2Index, -conductance);
                localTripletArray[index + 3] = Triplet(node2Index, node1Index, conductance);
            }
            else
            {
                nonLocalTripletArray.emplace_back(node2Index, node2Index, -conductance);
                nonLocalTripletArray.emplace_back(node2Index, node1Index, conductance);
            }
        }
    }
}
