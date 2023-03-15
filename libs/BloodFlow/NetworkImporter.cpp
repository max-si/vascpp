#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "BloodFlowEdgeNode.h"
#include "NetworkImporter.h"
#include "BloodFlowVessel.h"
#include "NetworkDescription.h"

// #include <vascular/core/BoundingBox.h>
// #include <vascular/core/Geometry.h>
// #include <vascular/core/GlobalDefs.h>

#include "VesselGenerator/coordinate.h"
#include "VesselGenerator/ImportExportCommon.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef vector<BloodFlowVessel> VesselVector;
typedef unsigned long ulong;
typedef Coordinate Point;

/**
 * \breif Imports and detrmines the number of vessels and starting position of partitions when multiple partitions are used
 *
 *
 * @param[in] fileId
 * @param[in, out] networkDescription
 * @param[in] numPartitions
 */
void MultiPartitionParallelImportPartitionInformation(hid_t fileId, NetworkDescription& networkDescription, int numPartitions)
{
    std::vector<long long> numVesselsPerPartition = ImportNumberOfVesselsPerPartitionData(fileId, numPartitions);
    //cout << "size of num VesselsPerPartition: " <<numVesselsPerPartition.size()<<"?=" << numPartitions<<endl;

    networkDescription.getNumberOfVesselsPerParition().assign(
        numVesselsPerPartition.begin(), numVesselsPerPartition.end());

    //calculate file start index for vessel import/export
    long long temp = 0;
    for (int i = 0; i < networkDescription.getLocalPartitionIndex(); i++)
    {
        temp += numVesselsPerPartition[i];
    }
    networkDescription.setGlobalVesselStartIndex(temp);

    //count & total number of vessels
    networkDescription.setTotalNumberOfVessels(std::accumulate(
        numVesselsPerPartition.begin(), numVesselsPerPartition.end(), 0LL));
}

/**
 * \breif Function that imports and sets parititon information for network
 *
 * Handles error reporting when number of partitions are equal to the of compute nodes
 * Also Handles the ability for a single compute node to run a partitioned file
 * @param[in] fileId
 * @param[in, out] networkDescription
 */
void ImportPartitionInformation(hid_t fileId, NetworkDescription& networkDescription)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    //cout << " IMport Partition Number: " << partitionNum << endl;

    int numPartitions = GetNumPartitionsInFile(fileId);

    //cout << numPartitions <<endl;
    if (mpiSize > 1 && numPartitions == mpiSize)
    {
        MultiPartitionParallelImportPartitionInformation(fileId, networkDescription, numPartitions);
    }
    else if (mpiSize == 1 && numPartitions == 1)
    {
        hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, GEOM_DATASET);
        hsize_t dims[2] = {0};
        hid_t dataspaceId = H5Dget_space(datasetId);
        H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
        networkDescription.setNumberOfVesselsInPartition(0, dims[0]);

        networkDescription.setGlobalVesselStartIndex(0);
        networkDescription.setTotalNumberOfVessels(dims[0]);

        H5Dclose(datasetId);
        H5Sclose(dataspaceId);
    }
    else if (mpiSize == 1 && numPartitions > 1)
    {
        std::vector<long long> numVesselsPerPartition = ImportNumberOfVesselsPerPartitionData(fileId, numPartitions);

        //calculate file start index for vessel import/export
        networkDescription.setGlobalVesselStartIndex(0);

        //count & total number of vessels
        networkDescription.setTotalNumberOfVessels(std::accumulate(
            numVesselsPerPartition.begin(), numVesselsPerPartition.end(), 0LL));
        networkDescription.setNumberOfVesselsInPartition(0, networkDescription.getTotalNumberOfVessels());
    }
    else
    {
#ifdef NDEBUG
        if (networkDescription.getNumberOfPartitions() != mpiSize)
        {
            cout << "Number of Partitions must be equal to the number of "
                    "processors. This might be relaxed in the future."
                 << endl;
            throw std::runtime_error(
                "Number of MPI Ranks Must Equal Number of Partitions");
        }
#endif
    }
    /*std::vector<long long>& numVesselsPerPartition = networkDescription.getNumberOfVesselsPerParition();
	numVesselsPerPartition.assign(numVesselsInPartitions, numVesselsInPartitions + numPartitions);*/

    //H5Sclose(memSpace);

    //delete[] numVesselsInPartitions;
}

/**
 * \breif Imports boundary conditions for nodes
 * @param[in] fileId
 * @param[in, out] networkDescription
 */
void ImportBcNodes(hid_t fileId, NetworkDescription& networkDescription)
{
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();
    hid_t datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    //cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    int numBoundaryConditions = dims[0];
    //import boundary condition Node IDs
    long long* nodeIndex = new long long[numBoundaryConditions];

    H5Dread(datasetId, H5T_NATIVE_LLONG, H5S_ALL, dataspace, dtplId, nodeIndex);

    H5Sclose(dataspace);
    H5Dclose(datasetId);
    // importing boundary condition type
    datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_TYPE_DATASET);
    dataspace = H5Dget_space(datasetId);

    short* bcType = new short[numBoundaryConditions];

    H5Dread(datasetId, H5T_NATIVE_SHORT, H5S_ALL, dataspace, dtplId, bcType);

    H5Sclose(dataspace);
    H5Dclose(datasetId);

    // import value of BC
    datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_VALUE_DATASET);
    dataspace = H5Dget_space(datasetId);

    double* bcValue = new double[numBoundaryConditions];

    H5Dread(datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace, dtplId, bcValue);

    H5Sclose(dataspace);
    H5Dclose(datasetId);

    //Store data

    for (int i = 0; i < numBoundaryConditions; i++)
    {
        short currentBcType = bcType[i];

        if (currentBcType < 0)
        {
            networkDescription.addBoundaryConditionNode(nodeIndex[i],
                NetworkDescription::BoundrayConditionType::Pressure,
                bcValue[i]);
        }
        if (currentBcType > 0)
        {
            networkDescription.addBoundaryConditionNode(nodeIndex[i],
                NetworkDescription::BoundrayConditionType::Flow, bcValue[i]);
        }
    }

    //delete all data
    H5Pclose(dtplId);
    delete[] nodeIndex;
    delete[] bcType;
    delete[] bcValue;
}

/**
 * \beief Imports Edge node indexes, partitions, and host index
 * @param[in] fileId
 * @param[in, out] networkDescription
 */
void ImportEdgeNodeInformation(
    hid_t fileId, NetworkDescription& networkDescription)
{
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();
    hid_t datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, EDGE_NODE_INDEX_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    //cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    int numEdgeNodes = dims[0];
    H5Sclose(dataspace);

    long long* index = new long long[numEdgeNodes];
    short* hostPartition = new short[numEdgeNodes];
    short* numAssociatedParts = new short[numEdgeNodes];

    H5Dread(datasetId, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dtplId, index);
    H5Dclose(datasetId);
    datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, EDGE_NODE_HOST_PARTITION_DATASET);
    H5Dread(datasetId, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, dtplId, hostPartition);
    H5Dclose(datasetId);

    datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, EDGE_NODE_NUM_ASSOCIATED_PARTITIONS_DATASET);
    H5Dread(datasetId, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, dtplId, numAssociatedParts);
    H5Dclose(datasetId);

    datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, EDGE_NODE_ASSOCIATED_PARTITIONS_DATASET);
    dataspace = H5Dget_space(datasetId);

    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    int associatedPartCount = dims[0];
    H5Sclose(dataspace);

    short* associatedParts = new short[associatedPartCount];

    H5Dread(datasetId, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, dtplId, associatedParts);
    H5Dclose(datasetId);

    //Add edge nodes to network description
    int associatedPartStartIndex = 0;

    for (int i = 0; i < numEdgeNodes; i++)
    {
        BloodFlowEdgeNode node(index[i]);
        node.setHostPartition(hostPartition[i]);
        for (int j = 0; j < numAssociatedParts[i]; j++)
        {
            node.addAssociatedPartition(associatedParts[associatedPartStartIndex + j]);
        }
        associatedPartStartIndex += numAssociatedParts[i];
        networkDescription.addEdgeNode(node);
    }

    delete[] index;
    delete[] hostPartition;
    delete[] numAssociatedParts;
    delete[] associatedParts;

    H5Pclose(dtplId);
}

/**
 * \breif Imports the non-vessel description of the vascular network
 * 
 * @param[in] fileId
 * @param[in, out] networkDesctiption
 */
void BuildNetworkDescription(
    hid_t fileId, NetworkDescription& networkDesctiption)
{
    //import partition information
    //cout << "IMport partition Information" << endl;
    ImportPartitionInformation(fileId, networkDesctiption);

    networkDesctiption.setNumLevels(RetrieveNumberOfLevelsFromFile(fileId));

    //import edge nodes
    //cout << "IMport Edge Nodes" << endl;
    if (networkDesctiption.getNumberOfPartitions() > 1)
        ImportEdgeNodeInformation(fileId, networkDesctiption);
    //import BC nodes
    //cout << "IMport BC Nodes" << endl;
    ImportBcNodes(fileId, networkDesctiption);
}

/**
 *  \breif Creates Vessel objects from imported data
 * @param[in] numVesselsToConstruct
 * @param[in, out] network
 * @param[in] geomData
 * @param[in] nodeData
 */
void ConstructVesselsFromData(long long numVesselsToConstruct,
    NetworkDescription& network, std::vector<double>& geomData,
    std::vector<long long>& nodeData)
{
    for (long long i = 0; i < numVesselsToConstruct; i++)
    {
        //cout << i << endl;
        vector<double>::size_type geomIndex = i * 7;
        vector<double>::size_type nodeIndex = i * 2;
        double length = CalculateLengthOfVessel(geomData[geomIndex],
            geomData[geomIndex + 1], geomData[geomIndex + 2],
            geomData[geomIndex + 3], geomData[geomIndex + 4],
            geomData[geomIndex + 5]);
        double radius = geomData[geomIndex + 6];
        long long node1 = nodeData[nodeIndex];
        long long node2 = nodeData[nodeIndex + 1];
        BloodFlowVessel currentVessel(length, radius, node1, node2);
        if (network.isBoundaryConditionNode(node1))
            currentVessel.setNode1IsBcNode(true);
        if (network.isBoundaryConditionNode(node2))
            currentVessel.setNode2IsBcNode(true);
        if (!network.isEdgeNode(node1))
        {
            currentVessel.setNode1IsLocal(true);
        }
        else
        {
            if (network.getLocalPartitionIndex() == network.getEdgeNode(node1).getHostPartition())
                currentVessel.setNode1IsLocal(true);
            else
                currentVessel.setNode1IsLocal(false);
        }
        if (!network.isEdgeNode(node2))
        {
            currentVessel.setNode2IsLocal(true);
        }
        else
        {
            if (network.getLocalPartitionIndex() == network.getEdgeNode(node2).getHostPartition())
                currentVessel.setNode2IsLocal(true);
            else
                currentVessel.setNode2IsLocal(false);
        }
        //cout << currentVessel << endl;
        network.addLocalVessel(currentVessel);
    }
}

/**
 * \breif Imports Vessel Data
 * 
 * Uses OpenMP to performed parallel block importing of vascular geometry data. 
 * This was necesary due to limitations in the HDF5 Library, it also prevents memory overruns of the deserializing process
 * 
 * @param[in] fileId
 * @param[in, out] networkDescription
 */
void ImportVesselData(hid_t fileId, NetworkDescription& networkDescription)
{
    networkDescription.clearLocalVesselVector();
    //cout << "calling vessel data importer" << endl;
    // open HDF datasets
    hid_t geomDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, GEOM_DATASET);
    hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);

    // determine import block size
    long long remainingNumberOfVesselsToExport = networkDescription.getLocalNumberOfVessels() % maxNumVesselsPerWrite;
    int numLoopsNeeded = networkDescription.getLocalNumberOfVessels() / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << remainingNumberOfVesselsToExport << " " << numLoopsNeeded << " " << globalNumLoopsNeeded <<" " << endl;
    //createRequiredArrays
    vector<double> geomImportArray;
    vector<double> geomBuildArray;

    vector<long long> nodeImportArray;
    vector<long long> nodeBuildArray;

    long long numVesselsToConstruct = 0;
    if (remainingNumberOfVesselsToExport > 0)
    {
        //cout << "Importing Remainder" << endl;
        ImportGeometryAndNodeDataBlock(geomDatasetId, nodeDatasetId,
            networkDescription.getGlobalVesselStartIndex(),
            remainingNumberOfVesselsToExport, geomBuildArray, nodeBuildArray);
        numVesselsToConstruct = remainingNumberOfVesselsToExport;
    }
    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
            //cout << "Importing Loop Block: " <<i << endl;
            long long vesselStartLocation = networkDescription.getGlobalVesselStartIndex() + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            {
                {
                    ImportGeometryAndNodeDataBlock(geomDatasetId, nodeDatasetId,
                        vesselStartLocation, maxNumVesselsPerWrite,
                        geomImportArray, nodeImportArray);
                }
                {
                    ConstructVesselsFromData(numVesselsToConstruct,
                        networkDescription, geomBuildArray, nodeBuildArray);
                }
            }

            std::swap(geomBuildArray, geomImportArray);
            std::swap(nodeBuildArray, nodeImportArray);
            numVesselsToConstruct = maxNumVesselsPerWrite;
        }
        else
        {
            //cout << "Null Importing Loop Block: " << i << endl;
            PerformNullImportOfGeometryAndNodes(geomDatasetId, nodeDatasetId);
        }
    }

    ConstructVesselsFromData(numVesselsToConstruct, networkDescription, geomBuildArray, nodeBuildArray);

    H5Dclose(geomDatasetId);
    H5Dclose(nodeDatasetId);
    //cout << "leaving vessel data importer" << endl;
}

/**
 * \breif Import utility to import damaged radii
 *
 * Needs to have the ability ot check if dataset exists added to the function
 * @param[in] fileId
 * @param[in, out] networkDescription
 */
void ImportDamagedRadii(hid_t fileId, NetworkDescription& networkDescription)
{
    // open HDF datasets
    hid_t datasetId = OpenHdfDataset(fileId, DAMAGE_GROUP_NAME, "/DAMAGED_RADII_ARRAY");

    vector<double> dataArray;

    CollectiveHdfBlockImport(datasetId,
        networkDescription.getGlobalVesselStartIndex(),
        networkDescription.getLocalNumberOfVessels(), dataArray);

    H5Dclose(datasetId);

    for (size_t i = 0; i < networkDescription.getLocalNumberOfVessels(); ++i)
    {
        networkDescription.getLocalVesselVector()[i].setDamagedRadius(dataArray[i]);
    }

    //cout << "leaving vessel data importer" << endl;
}

/**
 * \breif Determines matrix row bounds for local matrix partition
 * @param[in, out] network
 */
void ConstructRowBounds(NetworkDescription& network)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    vector<BloodFlowVessel>& vesselVector = network.getLocalVesselVector();
    long long minNodeIndex = std::numeric_limits<long long>::max();
    long long maxNodeIndex = std::numeric_limits<long long>::min();

    // Find min and max indexes of pure local junctions
    for (auto itr = vesselVector.begin(); itr != vesselVector.end(); ++itr)
    {
        if (itr->getNode1IsLocal())
        {
            if (minNodeIndex > itr->getNode1Index())
                minNodeIndex = itr->getNode1Index();
            if (maxNodeIndex < itr->getNode1Index())
                maxNodeIndex = itr->getNode1Index();
        }

        if (itr->getNode2IsLocal())
        {
            if (minNodeIndex > itr->getNode2Index())
                minNodeIndex = itr->getNode2Index();
            if (maxNodeIndex < itr->getNode2Index())
                maxNodeIndex = itr->getNode2Index();
        }
    }
    // Communicate all the min maxes to each other
    long long* recvBuffer = new long long[mpiSize * 2];
    long long sendBuffer[2];
    sendBuffer[0] = minNodeIndex;
    sendBuffer[1] = maxNodeIndex;

    MPI_Allgather(sendBuffer, 2, MPI_LONG_LONG_INT, recvBuffer, 2, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

    //add row pairs to network description
    for (int i = 0; i < mpiSize; ++i)
    {
        int index = i * 2;
        network.addRowBoundPairToVector(recvBuffer[index], recvBuffer[index + 1]);
    }

    delete[] recvBuffer;
}

/**
 * \breif Function that imports Network Description
 * @param[in] fileName
 * @param[in, out] networkDescription
 * @return HDF5 file handle
 */
hid_t BloodFlowNetworkImporter(
    std::string fileName, NetworkDescription& networkDescription)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    hid_t fileId = OpenHdfFile(fileName);

    networkDescription.setLocalPartitionIndex(mpiRank);
    //networkDescription.setLocalPartitionIndex(1);

    //build network information
    BuildNetworkDescription(fileId, networkDescription);
    //cout << networkDescription.getEdgeNodeSet().size() <<endl;
    if (networkDescription.getEdgeNodeSet().size() == 0 && networkDescription.getNumberOfPartitions() > 1)
    {
        cout << "no edge nodes in network that requires this information " << endl;
    }

    //import and construct vessels
    ImportVesselData(fileId, networkDescription);

    //! ImportDamagedRadii(fileId, networkDescription);

    //construct row bounds
    ConstructRowBounds(networkDescription);

    return fileId;
}