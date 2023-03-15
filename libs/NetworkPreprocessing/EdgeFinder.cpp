#include <algorithm>
#include <array>
#include <bitset>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

//#include <omp.h>
#include <mpi.h>

#include "EdgeFinder.h"
#include "EdgeFinderIO.h"
#include "EdgeFinderVessel.h"
#include "EdgeNode.h"
#include "PartitionInformation.h"

using std::array;
using std::bitset;
using std::cout;
using std::endl;
using std::map;
using std::set;
using std::string;
using std::vector;

//! I think this assumes all partitions have contiguous vessels
void MarkVesselLocality(long long vesselStartIndex, vector<EdgeFinderVessel>& vessels, map<long long, EdgeFinderVessel>& nonLocalVessels)
{
    long long lowerVesselIndexBound = vesselStartIndex;
    long long upperVesselIndexBound = vesselStartIndex + vessels.size();

    for (long long j = 0; j < vessels.size(); j++)
    {
        const vector<long long>& temp = vessels[j].getConnectedVessels();

        // std::cout << j << ": ";
        // for(int i = 0; i < temp.size(); i++) {
        //     std::cout << temp[i] << " ";
        // } std::cout << std::endl;
        
        bool allConnectedVesselsInRange = true;
        for (vector<long long>::const_iterator itr = temp.begin(); itr != temp.end(); itr++) {
            if (*itr < lowerVesselIndexBound || *itr >= upperVesselIndexBound) {
                allConnectedVesselsInRange = false;
                break;
            }
        }

        if (allConnectedVesselsInRange) {
            vessels[j].setVesselLocal(true);
        } else {
            vessels[j].setVesselLocal(false);
            nonLocalVessels.insert(std::make_pair(j + vesselStartIndex, vessels[j]));
        }
    }
}

std::unordered_map<long long, std::array<long long, 2>> CreateNonLocalVesselMap(
    const std::map<long long, EdgeFinderVessel>& nonLocalVessels,
    long long lowerIndexBound, long long upperIndexBound)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //Gather All vessel Data
    int numElementsToSend = nonLocalVessels.size() * 3;
    vector<int> numElementsOnEachRank(mpiSize);
    MPI_Allgather(&numElementsToSend, 1, MPI_INT, numElementsOnEachRank.data(),
        1, MPI_INT, MPI_COMM_WORLD);
    long long totalNumElementsToReceive = std::accumulate(
        numElementsOnEachRank.begin(), numElementsOnEachRank.end(), 0LL);

    vector<int> receiveOffsets = {0};

    std::partial_sum(numElementsOnEachRank.begin(),
        numElementsOnEachRank.end() - 1, std::back_inserter(receiveOffsets));
    vector<long long> dataToSend(numElementsToSend);
    vector<long long> receivedData(totalNumElementsToReceive);
    for (auto itr = nonLocalVessels.begin(); itr != nonLocalVessels.end();
         ++itr)
    {
        vector<long long>::size_type index =
            std::distance(nonLocalVessels.begin(), itr) * 3;
        dataToSend[index] = itr->first;
        dataToSend[index + 1] = itr->second.getStartNode();
        dataToSend[index + 2] = itr->second.getEndNode();
    }

    MPI_Allgatherv(dataToSend.data(), numElementsToSend, MPI_LONG_LONG_INT,
        receivedData.data(), numElementsOnEachRank.data(),
        receiveOffsets.data(), MPI_LONG_LONG_INT, MPI_COMM_WORLD);

    //Create Map
    std::unordered_map<long long, std::array<long long, 2>> vesselNodeMap;

    long long totalMapElements = totalNumElementsToReceive / 3 - numElementsToSend / 3;
    std::array<long long, 2> temp;
    for (vector<long long>::size_type i = 0; i < receivedData.size(); i += 3)
    {
        if (receivedData[i] < lowerIndexBound ||
            receivedData[i] >= upperIndexBound)
        {
            temp[0] = receivedData[i + 1];
            temp[1] = receivedData[i + 2];
            vesselNodeMap.insert(vesselNodeMap.end(), std::make_pair(receivedData[i], temp));
        }
    }

    return std::move(vesselNodeMap);
}

std::set<EdgeNode> CreateEdgeNodeSet(const vector<EdgeFinderVessel>& vessels,
    const map<long long, EdgeFinderVessel>& nonLocalVessels,
    const std::unordered_map<long long, std::array<long long, 2>> nonLocalVesselMap,
    const long long vesselStartIndex)
{
    std::set<EdgeNode> localEdgeNodes;
    long long lowerIndexBound = vesselStartIndex;
    long long upperIndexBound = vesselStartIndex + vessels.size();
    for (map<long long, EdgeFinderVessel>::const_iterator vesselItr = nonLocalVessels.begin(); vesselItr != nonLocalVessels.end(); ++vesselItr)
    {
        vector<EdgeNode> currentVesselNodes = vesselItr->second.getArrayOfNodes();

        vector<long long> currentVesselConnectedVessels = vesselItr->second.getConnectedVessels();

        for (vector<long long>::const_iterator connectedVesselItr = currentVesselConnectedVessels.begin();
             connectedVesselItr != currentVesselConnectedVessels.end(); ++connectedVesselItr)
        {
            if (*connectedVesselItr >= lowerIndexBound && *connectedVesselItr < upperIndexBound)
            {
                const EdgeFinderVessel& currentConnectedVessel = vessels[*connectedVesselItr - vesselStartIndex];
                long long matchingNode = currentConnectedVessel.FindMatchingNode(vesselItr->second);
                auto matchingEdgeNodeItr = std::find(currentVesselNodes.begin(), currentVesselNodes.end(), matchingNode);
                if (matchingEdgeNodeItr != currentVesselNodes.end())
                {
                    if (currentConnectedVessel.getVesselLocal())
                    {
                        currentVesselNodes.erase(matchingEdgeNodeItr);
                    }
                    else
                    {
                        matchingEdgeNodeItr->addLocalConnectedVessel(
                            *connectedVesselItr);
                    }
                }
            }
            else
            {
                const std::array<long long, 2>& temp = nonLocalVesselMap.find(*connectedVesselItr)->second;
                long long matchingNode = vesselItr->second.FindMatchingNode(temp);
                auto matchingEdgeNodeItr = std::find(currentVesselNodes.begin(), currentVesselNodes.end(), matchingNode);
                if (matchingEdgeNodeItr != currentVesselNodes.end())
                {
                    matchingEdgeNodeItr->addNonLocalConnectedVessel(*connectedVesselItr);
                }
            }
        }
        // add my index to remaining edge nodes
        for (auto& x : currentVesselNodes)
        {
            x.addLocalConnectedVessel(vesselItr->first);
            localEdgeNodes.insert(x);
        }
    }
    return std::move(localEdgeNodes);
}


//!!!!! Find edge node error
void FindLocalEdgeNodes(vector<EdgeFinderVessel>& vessels,
    long long vesselStartIndex, std::set<EdgeNode>& localEdgeNodes)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    map<long long, EdgeFinderVessel> nonLocalVessels;

    MarkVesselLocality(vesselStartIndex, vessels, nonLocalVessels);
    // cout << vessels.size() << endl;
    // for(int i = 0; i < nonLocalVessels.size(); i++){
    //     std::cout << nonLocalVessels[i] << std::endl;
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    if (!mpiRank) {cout << "Entering Vessel Node Map Creator" << endl;}
    std::unordered_map<long long, std::array<long long, 2>> nonLocalVesselMap =
        CreateNonLocalVesselMap(nonLocalVessels, vesselStartIndex, vesselStartIndex + vessels.size());
    //cout << "Leaving Vessel Node Map Creator" << endl;

    //if (!mpiRank) cout << "Create edge node Set" << endl; MPI_Barrier(MPI_COMM_WORLD);
    localEdgeNodes = CreateEdgeNodeSet(vessels, nonLocalVessels, nonLocalVesselMap, vesselStartIndex);
}

void AssignHostAndCommunicationPartitions(PartitionInformation& partInfo, std::set<EdgeNode>& localEdgeNodes)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long* startPointer = partInfo.numVesselsPerPartition;
    long long* endPointer = startPointer + partInfo.numPartitions + 1;

    for (auto& edgeNode : localEdgeNodes)
    {
        long long nodeNum = edgeNode.getNode();
        set<int> partitionNumVesselsMap;
        partitionNumVesselsMap.insert(partInfo.partitionNum);

        set<long long> nonLocalConnectedVessels = edgeNode.getNonLocalConnectedVessels();
        for (auto& connectedVessel : nonLocalConnectedVessels)
        {
            long long* upperBound = std::upper_bound(startPointer, endPointer, connectedVessel);
            int connectedVesselPartition = (upperBound - startPointer) - 1;
            if (!partitionNumVesselsMap.count(connectedVesselPartition))
            {
                partitionNumVesselsMap.insert(connectedVesselPartition);
            }
        }
        
        int hostPartition = *partitionNumVesselsMap.begin();
        int hostNumVessels = *partitionNumVesselsMap.begin();

        for (auto& partition : partitionNumVesselsMap)
        {
            const_cast<EdgeNode&>(edgeNode).addConnectedPartition(partition);
        }

        const_cast<EdgeNode&>(edgeNode).setHostPartition(hostPartition);    
			// for some reason i am getting a const iterator. No idea why so using hack to get around.
    }
}

void FindAndLabelEdgeNodes(hid_t fileId,
    std::vector<EdgeFinderVessel>& partitionVessels,
    PartitionInformation& partInfo, std::set<EdgeNode>& localEdgeNodes)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    partInfo.partitionNum = mpiRank;

    //import data
    MPI_Barrier(MPI_COMM_WORLD);
    ImportVesselsEdgeFinder(fileId, partInfo, partitionVessels);

    //find local edge nodes
    MPI_Barrier(MPI_COMM_WORLD);
    FindLocalEdgeNodes(partitionVessels, partInfo.startVesselIndex, localEdgeNodes);

    //determine communication pattern
    AssignHostAndCommunicationPartitions(partInfo, localEdgeNodes);
}

