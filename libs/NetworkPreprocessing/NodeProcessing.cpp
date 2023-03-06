#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

// #include <vascular/core/GlobalDefs.h>
#include "VesselGenerator/ImportExportCommon.h"
#include "VesselGenerator/ImportExportHdfBlock.h"
#include "EdgeFinder.h"
#include "EdgeFinderVessel.h"
#include "EdgeNode.h"
#include "NodeProcessing.h"
#include "PartitionInformation.h"

using std::cout;
using std::endl;

//! not used
void ExtractNodeIndexesAndCreateCommunicationMaps(
    std::vector<EdgeFinderVessel>& partitionVessels,
    std::vector<long long>& nodes,
    std::unordered_map<long long, long long>& newNodeIndexMap,
    std::set<EdgeNode>& localEdgeNodes, PartitionInformation& partInfo,
    std::map<int, std::set<long long>>& nonLocalNodeIndexesToSendMap,
    std::map<long long, long long>& nonLocalNodeIndexesToReceiveMap,
    long long indexOffset)
{
    long long newIndex = indexOffset;
    for (size_t i = 0; i < partitionVessels.size(); ++i)
    {
        std::vector<long long>::size_type arrayIndex = i * 2;
        long long node1 = partitionVessels[i].getStartNode();
        long long node2 = partitionVessels[i].getEndNode();
        nodes[arrayIndex] = node1;
        nodes[arrayIndex + 1] = node2;

        if (partitionVessels[i].getVesselLocal())
        {
            if (!newNodeIndexMap.count(node1))
            {
                newNodeIndexMap.insert(std::make_pair(node1, newIndex++));
            }
            if (!newNodeIndexMap.count(node2))
            {
                newNodeIndexMap.insert(std::make_pair(node2, newIndex++));
            }
        }
        else
        {
            if (localEdgeNodes.count(node1))
            {
                auto itr = localEdgeNodes.find(node1);
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    if (!newNodeIndexMap.count(node1))
                    {
                        newNodeIndexMap.insert(
                            std::make_pair(node1, newIndex++));
                    }
                }
                else
                {
                    nonLocalNodeIndexesToReceiveMap.insert(
                        std::make_pair(node1, itr->getHostPartition()));
                }
            }
            else
            {
                if (!newNodeIndexMap.count(node1))
                {
                    newNodeIndexMap.insert(std::make_pair(node1, newIndex++));
                }
            }

            if (localEdgeNodes.count(node2))
            {
                auto itr = localEdgeNodes.find(node2);
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    if (!newNodeIndexMap.count(node2))
                    {
                        newNodeIndexMap.insert(
                            std::make_pair(node2, newIndex++));
                    }
                }
                else
                {
                    nonLocalNodeIndexesToReceiveMap.insert(
                        std::make_pair(node2, itr->getHostPartition()));
                }
            }
            else
            {
                if (!newNodeIndexMap.count(node2))
                {
                    newNodeIndexMap.insert(std::make_pair(node2, newIndex++));
                }
            }
        }

        for (std::set<EdgeNode>::iterator itr = localEdgeNodes.begin();
             itr != localEdgeNodes.end(); itr++)
        {
            if (itr->getHostPartition() == partInfo.partitionNum)
            {
                const std::set<int>& associatedPartitions =
                    itr->getConnectedPartitions();
                for (std::set<int>::const_iterator partitionIndexItr =
                         associatedPartitions.begin();
                     partitionIndexItr != associatedPartitions.end();
                     ++partitionIndexItr)
                {
                    if (*partitionIndexItr != partInfo.partitionNum)
                    {
                        if (nonLocalNodeIndexesToSendMap.count(
                                *partitionIndexItr))
                        {
                            nonLocalNodeIndexesToSendMap
                                .find(*partitionIndexItr)
                                ->second.insert(itr->getNode());
                        }
                        else
                        {
                            std::set<long long> temp;
                            temp.insert(itr->getNode());
                            nonLocalNodeIndexesToSendMap.insert(
                                std::make_pair(*partitionIndexItr, temp));
                        }
                    }
                }
            }
        }
    }

    //int mpiSize, mpiRank;
    //MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    //MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //
    //for (int i = 0; i < mpiSize; i++)
    //{
    //	if (mpiRank == i)
    //	{
    //		/*cout << "Rank " << mpiRank << " nodes" << "\n";
    //		for (auto& x : nodes)
    //		{
    //
    //			cout << x << "\n";
    //		}
    //		cout << endl;*/
    //
    //		/*cout << "Rank " << mpiRank << " new Node map" << "\n";
    //		for (auto& x : newNodeIndexMap)
    //		{
    //
    //			cout << x.first << " "<< x.second << "\n";
    //		}
    //		cout << endl;*/
    //
    //		cout << "Rank " << mpiRank << " Non-Local Receive map" << "\n";
    //		for (auto& x : nonLocalNodeIndexesToReceiveMap)
    //		{
    //
    //			cout << x.first << " " << x.second << "\n";
    //		}
    //		cout << endl;
    //
    //		cout << "Rank " << mpiRank << " Non-Local Send map" << "\n";
    //		for (auto& x : nonLocalNodeIndexesToSendMap)
    //		{
    //
    //			cout << x.first << ": ";
    //			for(auto& y: x.second)
    //				cout << y <<" ";
    //
    //			cout << "\n";
    //		}
    //		cout << endl;
    //	}
    //	MPI_Barrier(MPI_COMM_WORLD);
    //}
}


void ExtractNodeIndexesAndCreateCommunicationMaps(
    std::vector<EdgeFinderVessel>& partitionVessels,
    std::vector<long long>& nodes,
    std::unordered_map<long long, long long>& newNodeIndexMap,
    std::set<EdgeNode>& localEdgeNodes, PartitionInformation& partInfo,
    std::vector<std::set<long long>>& nonLocalNodeIndexesToSendMap,
    std::vector<std::set<long long>>& nonLocalNodeIndexesToReceiveMap,
    long long indexOffset)
{
    long long newIndex = indexOffset;
    const long long numLoops = partitionVessels.size();
    for (size_t i = 0; i < numLoops; ++i)
    {
        std::vector<long long>::size_type arrayIndex = i * 2;
        const long long& node1 = partitionVessels[i].getStartNode();
        const long long& node2 = partitionVessels[i].getEndNode();
        nodes[arrayIndex] = node1;
        nodes[arrayIndex + 1] = node2;

        if (partitionVessels[i].getVesselLocal())
        {
            if (!newNodeIndexMap.count(node1))
            {
                newNodeIndexMap.insert(std::make_pair(node1, newIndex++));
            }
            if (!newNodeIndexMap.count(node2))
            {
                newNodeIndexMap.insert(std::make_pair(node2, newIndex++));
            }
        }
        else
        {
            if (localEdgeNodes.count(node1))
            {
                auto itr = localEdgeNodes.find(node1);
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    if (!newNodeIndexMap.count(node1))
                    {
                        newNodeIndexMap.insert(
                            std::make_pair(node1, newIndex++));
                    }
                }
                else
                {
                    nonLocalNodeIndexesToReceiveMap[itr->getHostPartition()]
                        .insert(node1);
                }
            }
            else
            {
                if (!newNodeIndexMap.count(node1))
                {
                    newNodeIndexMap.insert(std::make_pair(node1, newIndex++));
                }
            }

            if (localEdgeNodes.count(node2))
            {
                auto itr = localEdgeNodes.find(node2);
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    if (!newNodeIndexMap.count(node2))
                    {
                        newNodeIndexMap.insert(
                            std::make_pair(node2, newIndex++));
                    }
                }
                else
                {
                    nonLocalNodeIndexesToReceiveMap[itr->getHostPartition()]
                        .insert(node2);
                }
            }
            else
            {
                if (!newNodeIndexMap.count(node2))
                {
                    newNodeIndexMap.insert(std::make_pair(node2, newIndex++));
                }
            }
        }
    }

    for (std::set<EdgeNode>::iterator itr = localEdgeNodes.begin();
         itr != localEdgeNodes.end(); itr++)
    {
        if (itr->getHostPartition() == partInfo.partitionNum)
        {
            const std::set<int>& associatedPartitions =
                itr->getConnectedPartitions();
            for (std::set<int>::const_iterator partitionIndexItr =
                     associatedPartitions.begin();
                 partitionIndexItr != associatedPartitions.end();
                 ++partitionIndexItr)
            {
                if (*partitionIndexItr != partInfo.partitionNum)
                {
                    nonLocalNodeIndexesToSendMap[*partitionIndexItr].insert(
                        itr->getNode());
                }
            }
        }
    }
}

void RenumberMap(std::unordered_map<long long, long long>& newNodeIndexMap)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    long long localNumUniqueNodes = newNodeIndexMap.size();
    long long indexOffset = 0;
    MPI_Exscan(&localNumUniqueNodes, &indexOffset, 1, MPI_LONG_LONG_INT,
        MPI_SUM, MPI_COMM_WORLD);
    if (!mpiRank)
        indexOffset =
            0;    // sometimes mpi does not set this value for the host process

    for (std::unordered_map<long long, long long>::iterator itr =
             newNodeIndexMap.begin();
         itr != newNodeIndexMap.end(); itr++)
    {
        long long index = std::distance(newNodeIndexMap.begin(), itr);
        itr->second = index + indexOffset;
    }

    /*for (int i = 0; i < mpiSize; i++)
	{
		if (mpiRank == i)
		{
			cout << localNumUniqueNodes << " " << indexOffset << endl;
			cout << "Rank " << mpiRank << " new Node map" << "\n";
			for (auto& x : newNodeIndexMap)
			{
				cout << x.first << " " << x.second << "\n";
			}
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
}

void CommunicateEdgeNodeMaps(
    std::map<int, std::set<long long>>& nonLocalNodeIndexesToSendMap,
    std::unordered_map<long long, long long>& newNodeIndexMap,
    std::map<long long, long long>& nonLocalNodeIndexesToReceiveMap)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    std::vector<int> sendCount(mpiSize);
    std::vector<long long> sendData;
    for (auto itr = nonLocalNodeIndexesToSendMap.begin();
         itr != nonLocalNodeIndexesToSendMap.end(); ++itr)
    {
        for (auto& nodeNum : itr->second)
        {
            sendCount[itr->first]++;
            long long newNodeIndex = newNodeIndexMap.find(nodeNum)->second;
            sendData.push_back(newNodeIndex);
        }
    }

    std::vector<int> sendOffsets(mpiSize);
    std::partial_sum(
        sendCount.begin(), sendCount.end() - 1, sendOffsets.begin() + 1);

    std::vector<int> receiveCount(mpiSize);
    std::vector<int> receiveOffsets(mpiSize);

    int totalToReceive = 0;
    for (auto itr = nonLocalNodeIndexesToReceiveMap.begin();
         itr != nonLocalNodeIndexesToReceiveMap.end(); itr++)
    {
        receiveCount[itr->second]++;
        ++totalToReceive;
    }
    std::partial_sum(receiveCount.begin(), receiveCount.end() - 1,
        receiveOffsets.begin() + 1);
    std::vector<long long> receiveData(totalToReceive);

    MPI_Alltoallv(sendData.data(), sendCount.data(), sendOffsets.data(),
        MPI_LONG_LONG_INT, receiveData.data(), receiveCount.data(),
        receiveOffsets.data(), MPI_LONG_LONG_INT, MPI_COMM_WORLD);

    //for (int i = 0; i < mpiSize; i++)
    {
        //if (mpiRank == i)
        {
            //cout << "Rank " << i << " Mapping Received Data to Node Indexes" << endl;
            for (auto itr = nonLocalNodeIndexesToReceiveMap.begin();
                 itr != nonLocalNodeIndexesToReceiveMap.end(); itr++)
            {
                long long index =
                    std::distance(nonLocalNodeIndexesToReceiveMap.begin(), itr);
                //cout << "index: " << index << " Orig Number: " << itr->first << " New Value "<< receiveData[index] << endl;
                newNodeIndexMap.insert(
                    std::make_pair(itr->first, receiveData[index]));
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }

    /*for (int i = 0; i < mpiSize; i++)
	{
		if (mpiRank == i)
		{
			cout << "Rank " << mpiRank << " Node map\n\t";
			for (auto& x : newNodeIndexMap)
			{
				cout << x.first << " -> " << x.second << "\n\t";
			}
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	for (int i = 0; i < mpiSize; i++)
	{
		if (mpiRank == i)
		{
			cout << "Rank " << mpiRank << " send Counts\n";
			for (auto& x : sendCount)
			{
				cout << x << " ";
			}
			cout <<"\n";
			cout << "Rank " << mpiRank << " send Offsets\n";
			for (auto& x : sendOffsets)
			{
				cout << x << " ";
			}
			cout << "\n";
			cout << "Rank " << mpiRank << " send data\n";
			for (auto& x : sendData)
			{
				cout << x << " ";
			}
			cout << "\n";
			cout << "Rank " << mpiRank << " receive Counts\n";
			for (auto& x : receiveCount)
			{
				cout << x << " ";
			}
			cout << "\n";
			cout << "Rank " << mpiRank << " receive Offsets\n";
			for (auto& x : receiveOffsets)
			{
				cout << x << " ";
			}
			cout << "\n";
			cout << "Rank " << mpiRank << " receive data\n";
			for (auto& x : receiveData)
			{
				cout << x << " ";
			}
			cout << "\n";
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
}

void CommunicateEdgeNodeMaps(
    std::vector<std::set<long long>>& nonLocalNodeIndexesToSendMap,
    std::unordered_map<long long, long long>& newNodeIndexMap,
    std::vector<std::set<long long>>& nonLocalNodeIndexesToReceiveMap)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    std::vector<int> sendCount(mpiSize);
    std::vector<long long> sendData;

    for (int i = 0; i < nonLocalNodeIndexesToSendMap.size(); ++i)
    {
        sendCount[i] += nonLocalNodeIndexesToSendMap[i].size();
        for (auto& nodeNum : nonLocalNodeIndexesToSendMap[i])
        {
            long long newNodeIndex = newNodeIndexMap.find(nodeNum)->second;
            sendData.push_back(newNodeIndex);
        }
    }

    std::vector<int> sendOffsets(mpiSize);
    std::partial_sum(
        sendCount.begin(), sendCount.end() - 1, sendOffsets.begin() + 1);

    std::vector<int> receiveCount(mpiSize);
    std::vector<int> receiveOffsets(mpiSize);

    long long totalToReceive = 0;
    for (int i = 0; i < nonLocalNodeIndexesToReceiveMap.size(); ++i)
    {
        receiveCount[i] = nonLocalNodeIndexesToReceiveMap[i].size();
        totalToReceive += nonLocalNodeIndexesToReceiveMap[i].size();
    }
    std::partial_sum(receiveCount.begin(), receiveCount.end() - 1,
        receiveOffsets.begin() + 1);
    std::vector<long long> receiveData(totalToReceive);

    MPI_Alltoallv(sendData.data(), sendCount.data(), sendOffsets.data(),
        MPI_LONG_LONG_INT, receiveData.data(), receiveCount.data(),
        receiveOffsets.data(), MPI_LONG_LONG_INT, MPI_COMM_WORLD);

    long long dataIndex = 0;
    for (int i = 0; i < nonLocalNodeIndexesToReceiveMap.size(); ++i)
    {
        for (auto itr = nonLocalNodeIndexesToReceiveMap[i].begin(); itr != nonLocalNodeIndexesToReceiveMap[i].end(); ++itr) {
            newNodeIndexMap.insert(std::make_pair(*itr, receiveData[dataIndex++]));
		}
    }
}

void ApplyMapToNodes(std::vector<long long>& nodes, std::unordered_map<long long, long long>& map)
{
//#pragma omp parallel for
    for (std::vector<long long>::iterator nodeItr = nodes.begin(); nodeItr != nodes.end(); ++nodeItr)
    {
        *nodeItr = map.find(*nodeItr)->second;
    }
}

void ApplyMapToEdgeNodes(std::set<EdgeNode>& nodes, std::unordered_map<long long, long long>& map)
{
    std::set<EdgeNode> tempSet;
    for (std::set<EdgeNode>::iterator nodeItr = nodes.begin(); nodeItr != nodes.end(); ++nodeItr)
    {
        EdgeNode tempNode = *nodeItr;
        long long newIndex = map.find(tempNode.getNode())->second;
        tempNode.setNode(newIndex);
        tempSet.insert(tempSet.end(), tempNode);
    }
    nodes.swap(tempSet);
}

std::vector<long long> CollectiveBcNodeImport(hid_t fileId)
{
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();
    hid_t datasetId = OpenHdfDataset(
        fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    hsize_t numDims = H5Sget_simple_extent_ndims(dataspace);
    //cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    int numBoundaryConditions = dims[0];
    //import boundary condition Node IDs

    std::vector<long long> bcNodeIndexes(numBoundaryConditions);

    H5Dread(datasetId, H5T_NATIVE_LLONG, H5S_ALL, dataspace, dtplId,
        bcNodeIndexes.data());

    H5Sclose(dataspace);
    H5Dclose(datasetId);

    //delete all data
    H5Pclose(dtplId);
    return bcNodeIndexes;
}

void CollectiveBcNodeExport(hid_t fileId, std::vector<long long>& outputData, std::vector<int>& offsets)
{
    if (outputData.size() == 0)
        return;
    hid_t dtplId = CreateIndependentDataTransferPropertiesList();
    hid_t datasetId = OpenHdfDataset(
        fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    hsize_t numDims = H5Sget_simple_extent_ndims(dataspace);
    //cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    int numBoundaryConditions = dims[0];
    dims[0] = outputData.size();
    hid_t memspaceId = H5Screate_simple(1, dims, NULL);
    //import boundary condition Node IDs

    hsize_t start, block = 1, stride = 1;
    dims[0] = 1;
    start = offsets[0];

    //do first element to select hyperspace
    H5Sselect_hyperslab(
        dataspace, H5S_SELECT_SET, &start, &stride, dims, &block);

    for (int i = 1; i < outputData.size(); ++i)
    {
        start = offsets[i];
        H5Sselect_hyperslab(
            dataspace, H5S_SELECT_OR, &start, &stride, dims, &block);
    }

    H5Dwrite(datasetId, H5T_NATIVE_LLONG, memspaceId, dataspace, dtplId,
        outputData.data());

    H5Sclose(dataspace);
    H5Sclose(memspaceId);
    H5Dclose(datasetId);

    //delete all data
    H5Pclose(dtplId);
}

void ApplyMapToBcNodes(std::vector<long long>& nodes,
    std::unordered_map<long long, long long>& map,
    std::vector<long long>& outputData, std::vector<int>& offsets)
{
    for (int i = 0; i < nodes.size(); ++i)
    {
        if (map.count(nodes[i]))
        {
            offsets.push_back(i);
            long long newIndex = map[nodes[i]];
            outputData.push_back(newIndex);
        }
    }
}

void ProcessBcNodes(hid_t fileId, std::unordered_map<long long, long long>& map)
{
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<long long> bcNodes = ImportAndBcastSmallDataset<long long>(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET, MPI_LONG_LONG_INT);

    std::vector<int> newNodeOffsets;
    std::vector<long long> newNodeIndexes;

    MPI_Barrier(MPI_COMM_WORLD);
    ApplyMapToBcNodes(bcNodes, map, newNodeIndexes, newNodeOffsets);
    
    MPI_Barrier(MPI_COMM_WORLD);
    CollectiveBcNodeExport(fileId, newNodeIndexes, newNodeOffsets);
}

void ExportNodeData(hid_t fileId, PartitionInformation& partInfo, std::vector<long long>& nodeData)
{
    //cout << "calling vessel data importer" << endl;
    // open HDF datasets
    hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
    // determine import block size
    long long remainingNumberOfVesselsToExport =
        partInfo.numVesselsInPartition % maxNumVesselsPerWrite;
    int numLoopsNeeded = partInfo.numVesselsInPartition / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(
        remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << globalNumLoopsNeeded << " " << numLoopsNeeded <<" "<< remainingNumberOfVesselsToExport <<endl;
    //createRequiredArrays

    std::vector<long long> nodeExportArray;
    std::vector<long long> nodeDataArray;

    long long numVesselsToExport = 0;
    if (remainingNumberOfVesselsToExport > 0)
    {
        nodeExportArray.assign(nodeData.begin(),
            nodeData.begin() + (remainingNumberOfVesselsToExport * 2));
        numVesselsToExport = remainingNumberOfVesselsToExport;
    }
    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
//#pragma omp parallel sections num_threads(2)
            {
//#pragma omp section
                {
                    long long vesselStartLocation;
                    if (i == 0)
                    {
                        vesselStartLocation = partInfo.startVesselIndex;
                    }
                    else
                        vesselStartLocation = partInfo.startVesselIndex +
                            remainingNumberOfVesselsToExport +
                            (i - 1) * maxNumVesselsPerWrite;

                    //cout << "Block Export Loop " << i <<" " << vesselStartLocation << " " << numVesselsToExport << endl;
                    CollectiveHdfBlockExport<long long, 2, 2>(nodeDatasetId,
                        vesselStartLocation, numVesselsToExport,
                        nodeExportArray);
                }
//#pragma omp section
                {
                    long long start = (remainingNumberOfVesselsToExport * 2) +
                        i * (maxNumVesselsPerWrite * 2);
                    nodeDataArray.assign(nodeData.begin() + start,
                        nodeData.begin() + start + (maxNumVesselsPerWrite * 2));
                }
            }

            std::swap(nodeDataArray, nodeExportArray);
            numVesselsToExport = maxNumVesselsPerWrite;
        }
        else
        {
            CollectiveHdfNullExport<long long, 2, 2>(nodeDatasetId);
        }
    }
    long long vesselStartLocation;
    if (numLoopsNeeded == 0)
    {
        vesselStartLocation = partInfo.startVesselIndex;
        //cout << "Block Export No Loops " << vesselStartLocation << " " << numVesselsToExport << endl;
        CollectiveHdfBlockExport<long long, 2, 2>(nodeDatasetId,
            vesselStartLocation, numVesselsToExport, nodeExportArray);
    }
    else
    {
        vesselStartLocation = partInfo.startVesselIndex +
            remainingNumberOfVesselsToExport +
            (numLoopsNeeded - 1) * maxNumVesselsPerWrite;
        //cout << "Block Export After Loops " << vesselStartLocation << " " << numVesselsToExport << endl;
        CollectiveHdfBlockExport<long long, 2, 2>(nodeDatasetId,
            vesselStartLocation, numVesselsToExport, nodeExportArray);
    }

    H5Dclose(nodeDatasetId);
    //cout << "leaving vessel data importer" << endl;
}

long long DetermineNewNodeIndexOffsets(
    const std::vector<EdgeFinderVessel>& partitionVessels,
    const std::set<EdgeNode>& localEdgeNodes,
    const PartitionInformation& partInfo)
{
    std::unordered_set<long long> nodeSet;
    nodeSet.reserve(partitionVessels.size());
    //std::set<long long> nodeSet;
    long long numLoops = partitionVessels.size();
    const auto endItr = localEdgeNodes.end();
    //double totalTime = -omp_get_wtime(), insertTime = 0;
    for (int i = 0; i < numLoops; ++i)
    {
        const long long& node1 = partitionVessels[i].getStartNode();
        const long long& node2 = partitionVessels[i].getEndNode();
        if (!nodeSet.count(node1))
        {
            const auto itr = localEdgeNodes.find(node1);
            if (itr != endItr)
            {
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    //insertTime -= omp_get_wtime();
                    nodeSet.insert(node1);
                    //insertTime += omp_get_wtime();
                }
            }
            else
            {
                //insertTime -= omp_get_wtime();
                nodeSet.insert(node1);
                //insertTime += omp_get_wtime();
            }
        }

        if (!nodeSet.count(node2))
        {
            const auto itr = localEdgeNodes.find(node2);
            if (itr != endItr)
            {
                if (itr->getHostPartition() == partInfo.partitionNum)
                {
                    //insertTime -= omp_get_wtime();
                    nodeSet.insert(node2);
                    //insertTime += omp_get_wtime();
                }
            }
            else
            {
                //insertTime -= omp_get_wtime();
                nodeSet.insert(node2);
                //insertTime += omp_get_wtime();
            }
        }
    }
    //totalTime += omp_get_wtime();
    //cout << totalTime << " " << insertTime << " " << totalTime - insertTime << endl;
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long localNumUniqueNodes = nodeSet.size();
    long long indexOffset = 0;
    MPI_Exscan(&localNumUniqueNodes, &indexOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    if (!mpiRank) {indexOffset = 0;}
    return indexOffset;
}

long long DetermineNewNodeIndexOffsets2(
    const std::vector<EdgeFinderVessel>& partitionVessels,
    const std::set<EdgeNode>& localEdgeNodes,
    const PartitionInformation& partInfo)
{
    std::unordered_set<long long> nodeSet;
    nodeSet.reserve(partitionVessels.size());
    std::unordered_set<long long> nonLocalEdgeNodes;
    nonLocalEdgeNodes.reserve(localEdgeNodes.size());

    for (auto& x : localEdgeNodes)
    {
        if (x.getHostPartition() != partInfo.partitionNum)
            nonLocalEdgeNodes.insert(x.getNode());
    }

    long long numLoops = partitionVessels.size();
    const auto endItr = localEdgeNodes.end();
    
    for (int i = 0; i < numLoops; ++i)
    {
        const long long& node1 = partitionVessels[i].getStartNode();
        const long long& node2 = partitionVessels[i].getEndNode();
        if (!nodeSet.count(node1) && !nonLocalEdgeNodes.count(node1))
        {
            nodeSet.insert(node1);
        }

        if (!nodeSet.count(node2) && !nonLocalEdgeNodes.count(node2))
        {
            nodeSet.insert(node2);
        }
    }

    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long localNumUniqueNodes = nodeSet.size();
    long long indexOffset = 0;
    MPI_Exscan(&localNumUniqueNodes, &indexOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    if (!mpiRank) {indexOffset = 0;}
    return indexOffset;
}

void NodeProcessing(hid_t fileId) {
	//std::cout << "Entering Node Processing" << std::endl;
	int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    std::vector<EdgeFinderVessel> partitionVessels;
    PartitionInformation partInfo;
    std::set<EdgeNode> localEdgeNodes;

	FindAndLabelEdgeNodes(fileId, partitionVessels, partInfo, localEdgeNodes);

	MPI_Barrier(MPI_COMM_WORLD);
    
    std::vector<long long> nodes(2 * partitionVessels.size());
    std::unordered_map<long long, long long> newNodeIndexMap;
    std::vector<std::set<long long>> nonLocalNodeIndexesToSendMap(partInfo.numPartitions);
    std::vector<std::set<long long>> nonLocalNodeIndexesToReceiveMap(partInfo.numPartitions);
    
    long long indexOffset = DetermineNewNodeIndexOffsets2(partitionVessels, localEdgeNodes, partInfo);

	ExtractNodeIndexesAndCreateCommunicationMaps(partitionVessels, nodes,
        newNodeIndexMap, localEdgeNodes, partInfo, nonLocalNodeIndexesToSendMap,
        nonLocalNodeIndexesToReceiveMap, indexOffset);

	partitionVessels.clear();

	CommunicateEdgeNodeMaps(nonLocalNodeIndexesToSendMap, newNodeIndexMap, nonLocalNodeIndexesToReceiveMap);
    
    nonLocalNodeIndexesToReceiveMap.clear();
    nonLocalNodeIndexesToSendMap.clear();
    
    ApplyMapToEdgeNodes(localEdgeNodes, newNodeIndexMap);
    ApplyMapToNodes(nodes, newNodeIndexMap);

	ProcessBcNodes(fileId, newNodeIndexMap);

	newNodeIndexMap.clear();

	MPI_Barrier(MPI_COMM_WORLD);
    
    ExportNodeData(fileId, partInfo, nodes);
    ExportEdgeNodesToFile(fileId, localEdgeNodes, partInfo);
}