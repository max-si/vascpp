#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "VesselGenerator/ImportExportCommon.h"
#include "core/Timsort.h"
#include "EdgeFinderIO.h"
#include "EdgeFinderVessel.h"
#include "EdgeNode.h"
#include "PartitionInformation.h"

using gfx::timsort;
using std::cout;
using std::endl;
using std::set;
using std::string;
using std::vector;

void ImportNumVesselsInPartition(hid_t fileId, int partitionNum, PartitionInformation& partInfo)
{
    //cout << " IMport Partition Number: " << partitionNum << endl;
    hid_t datasetId = OpenHdfDataset(fileId, "/", NUM_VESSELS_PER_PARTITION_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    hsize_t numDims = H5Sget_simple_extent_ndims(dataspace);
    //cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    partInfo.numPartitions = dims[0];
    long long* numVesselsInPartitions = new long long[dims[0] + 1];
    hid_t memSpace = H5Screate_simple(1, dims, NULL);
    numVesselsInPartitions[0] = 0;

    H5Dread(datasetId, H5T_NATIVE_LLONG, memSpace, dataspace, H5P_DEFAULT, numVesselsInPartitions + 1);
    partInfo.numVesselsInPartition = numVesselsInPartitions[partitionNum + 1];

    long long temp = 0;
    for (int i = 0; i < partitionNum; i++)
        temp += numVesselsInPartitions[i + 1];
    partInfo.startVesselIndex = temp;

    for (int i = 1; i < dims[0] + 1; i++)

        numVesselsInPartitions[i] += numVesselsInPartitions[i - 1];

    partInfo.numVesselsPerPartition = numVesselsInPartitions;

    H5Sclose(memSpace);
    H5Sclose(dataspace);
    H5Dclose(datasetId);
    //delete[] numVesselsInPartitions;
}

void ImportConnectedVesselIndexes(hid_t fileId, int partitionNum, PartitionInformation& partInfo)  
{
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    herr_t status;
    hid_t dataTransferPropertyList = CreateCollectiveDataTransferPropertiesList();

    //create hyperslabs in dataset and file

    //	cout << "MPI_rank" << mpiRank << " reporting from File import" << endl;

    hsize_t dims[2], count[2], offset[2];
    long long numVesselsOnNode = partInfo.numVesselsInPartition;

    count[0] = numVesselsOnNode;
    count[1] = 1;
    hid_t memSpace = H5Screate_simple(2, count, NULL);

    offset[0] = partInfo.startVesselIndex;
    offset[1] = 0;

    //import number of connected vessels for node
    hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NUM_CONNECTED_VESSELS_DATASET);

    short* numConnectedVessels = new short[numVesselsOnNode];
    partInfo.numConnectedVesselsArray = numConnectedVessels;

    hid_t fileSpace = H5Dget_space(datasetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dread(datasetId, H5T_NATIVE_SHORT, memSpace, fileSpace, dataTransferPropertyList, numConnectedVessels);
    H5Sclose(fileSpace);
    H5Sclose(memSpace);
    H5Dclose(datasetId);

    //cout << "Import Connected Vessels (after import num connected Vessels): " << numConnectedVessels[0] << endl;

    //calculate total elements needed for import of connected vessels
    long long nodeNumConnectedVessels = 0;
    for (long long i = 0; i < numVesselsOnNode; i++)
    {
        nodeNumConnectedVessels += numConnectedVessels[i];
    }
    partInfo.numConnectedVessels = nodeNumConnectedVessels;
    //MPI_Ex_scan to get start position of import
    long long startPosition = 0;
    if (mpiSize > 1)
        MPI_Exscan(&nodeNumConnectedVessels, &startPosition, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpiRank == 0)
        startPosition = 0;

    partInfo.connectedVesselStartIndex = startPosition;
}

void ImportDataAndCreateVessels(PartitionInformation& partInfo, hid_t fileId, int mpiRank, std::vector<EdgeFinderVessel>& vessels) 
{
    long long remainingNumberOfVesselsToExport = partInfo.numVesselsInPartition % maxNumVesselsPerWrite;
    int numLoopsNeeded = partInfo.numVesselsInPartition / maxNumVesselsPerWrite;
    bool localNumRemainingVessels = remainingNumberOfVesselsToExport;
    bool anyProcessHasRemainingVessels = false;
    MPI_Allreduce(&localNumRemainingVessels, &anyProcessHasRemainingVessels, 1, MPI_CHAR, MPI_LOR, MPI_COMM_WORLD);
    if (anyProcessHasRemainingVessels && !localNumRemainingVessels)
    {
        remainingNumberOfVesselsToExport = maxNumVesselsPerWrite;
        numLoopsNeeded -= 1;
    }
    int globalNumLoopsNeeded = 0;
    MPI_Allreduce(&numLoopsNeeded, &globalNumLoopsNeeded, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    long long* nodeArray = new long long[maxNumVesselsPerWrite * 2];
    long long* connectedVesselIndexArray = static_cast<long long*>(malloc(sizeof(long long) * maxNumVesselsPerWrite * 5));

    hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
    hid_t connectedVesselDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);

    hid_t nodeMemSpaceId, nodeDataSpaceId;
    hsize_t nodeDims[2] = {0};
    hsize_t nodeStart[2] = {0}, nodeCount[2] = {0}, nodeBlockStride[2] = {0};
    hid_t connectedVesselMemSpaceId, connectedVesselDataSpaceId;
    hsize_t connectedVesselDims[2] = {0};
    hsize_t connectedVesselStart[2] = {0}, connectedVesselCount[2] = {0}, connectedVesselBlockStride[2] = {0};

    nodeBlockStride[0] = 1;
    nodeBlockStride[1] = 1;

    connectedVesselBlockStride[0] = 1;
    connectedVesselBlockStride[1] = 1;

    nodeDataSpaceId = H5Dget_space(nodeDatasetId);
    connectedVesselDataSpaceId = H5Dget_space(connectedVesselDatasetId);

    //hid_t ftplId = CreateIndependentDataTransferPropertiesList();
    hid_t ftplId = CreateCollectiveDataTransferPropertiesList();

    //if (!mpiRank) cout << "Import Odd Block" << endl;
    if (remainingNumberOfVesselsToExport > 0)
    {
        //define nodespaces
        nodeCount[0] = remainingNumberOfVesselsToExport;
        nodeCount[1] = 2;
        nodeMemSpaceId = H5Screate_simple(2, nodeCount, NULL);

        nodeStart[0] = partInfo.startVesselIndex;
        nodeStart[1] = 0;
        //cout << "remainingNumberOfVesselsToExport: " << remainingNumberOfVesselsToExport << endl;

        H5Sselect_hyperslab(nodeDataSpaceId, H5S_SELECT_SET, nodeStart, nodeBlockStride, nodeCount, nodeBlockStride);

        //import nodes
        H5Dread(nodeDatasetId, H5T_NATIVE_LLONG, nodeMemSpaceId, nodeDataSpaceId, ftplId, nodeArray);
        //H5Sclose(nodeDataSpaceId);

        /*cout << endl;
		for (long long i = 0; i < remainingNumberOfVesselsToExport*2; i++)
		cout << nodeArray[i] << endl;
		cout << endl;*/

        //construct base vessels
        for (long long i = 0; i < remainingNumberOfVesselsToExport; i++)
        {
            long long nodeIndex = i * 2;
            long long startNode = nodeArray[nodeIndex];
            long long endNode = nodeArray[nodeIndex + 1];
            short partitionNumber = mpiRank;
            //cout << i << " " << nodeIndex << " " <<startNode << " " <<endNode <<endl;
            vessels.emplace_back(startNode, endNode, partitionNumber);
        }

        //possibly parallel task with construct base vessels
        //define connected Vessel spaces
        connectedVesselCount[0] = 0;
        connectedVesselCount[1] = 1;
        for (long long i = 0; i < remainingNumberOfVesselsToExport; i++)
        {
            //cout << partInfo.numConnectedVesselsArray[i] << endl;
            connectedVesselCount[0] += partInfo.numConnectedVesselsArray[i];
        }
        //cout << "connected Vessel count " << connectedVesselCount[0] << endl;
        connectedVesselMemSpaceId = H5Screate_simple(2, connectedVesselCount, NULL);

        connectedVesselStart[0] = partInfo.connectedVesselStartIndex;
        connectedVesselStart[1] = 0;

        H5Sselect_hyperslab(connectedVesselDataSpaceId, H5S_SELECT_SET,
            connectedVesselStart, connectedVesselBlockStride,
            connectedVesselCount, connectedVesselBlockStride);
        // import connected vessel data
        H5Dread(connectedVesselDatasetId, H5T_NATIVE_LLONG,
            connectedVesselMemSpaceId, connectedVesselDataSpaceId, ftplId,
            connectedVesselIndexArray);
        //H5Sclose(connectedVesselDataSpaceId);
        // add connected vessel data to vessels
        long long* startPointer = connectedVesselIndexArray;
        for (long long i = 0; i < remainingNumberOfVesselsToExport; i++)
        {
            //cout << "hi" << endl;
            long long* endPointer = startPointer + partInfo.numConnectedVesselsArray[i];

            vector<long long> temp(startPointer, endPointer);
            //std::sort(temp.begin(), temp.end());
            timsort(temp.begin(), temp.end());

            vessels[i].setConnectedVessels(temp);
            startPointer = endPointer;
        }
    }

    // import rest of vessels in loop

    long long connectedVesselStartOffset = partInfo.connectedVesselStartIndex + connectedVesselCount[0];
    //cout << "Connected Vessel Parameters " << partInfo.numConnectedVessels << " " << partInfo.connectedVesselStartIndex << " " << connectedVesselStartOffset << endl;
    //if (!mpiRank) cout << "import loop blocks" << endl;
    for (int k = 0; k < globalNumLoopsNeeded; k++)
    {
        if (k < numLoopsNeeded)
        {
            //define nodespaces
            nodeCount[0] = maxNumVesselsPerWrite;
            nodeCount[1] = 2;
            nodeMemSpaceId = H5Screate_simple(2, nodeCount, NULL);

            nodeStart[0] = partInfo.startVesselIndex + remainingNumberOfVesselsToExport + k * maxNumVesselsPerWrite;
            nodeStart[1] = 0;

            H5Sselect_hyperslab(nodeDataSpaceId, H5S_SELECT_SET, nodeStart, nodeBlockStride, nodeCount, nodeBlockStride);

            //import nodes
            H5Dread(nodeDatasetId, H5T_NATIVE_LLONG, nodeMemSpaceId, nodeDataSpaceId, ftplId, nodeArray);
            H5Sclose(nodeMemSpaceId);
            //construct base vessels
            for (long long i = 0; i < maxNumVesselsPerWrite; i++)
            {
                long long nodeIndex = i * 2;
                long long startNode = nodeArray[nodeIndex];
                long long endNode = nodeArray[nodeIndex + 1];
                short partitionNumber = mpiRank;
                //cout << i << " " << nodeIndex << " " << startNode << " " << endNode << endl;
                vessels.emplace_back(startNode, endNode, partitionNumber);
            }

            //possibly parallel task with construct base vessels
            //define connected Vessel spaces
            connectedVesselCount[0] = 0;
            connectedVesselCount[1] = 1;
            for (long long i = 0; i < maxNumVesselsPerWrite; i++)
            {
                long long vesselIndex = remainingNumberOfVesselsToExport + k * maxNumVesselsPerWrite + i;
                connectedVesselCount[0] += partInfo.numConnectedVesselsArray[vesselIndex];
            }
            connectedVesselMemSpaceId = H5Screate_simple(2, connectedVesselCount, NULL);

            connectedVesselStart[0] = connectedVesselStartOffset;
            connectedVesselStart[1] = 0;
            //cout << k << " " << connectedVesselStart[0] << " " << connectedVesselCount[0] << endl;
            H5Sselect_hyperslab(connectedVesselDataSpaceId, H5S_SELECT_SET,
                connectedVesselStart, connectedVesselBlockStride,
                connectedVesselCount, connectedVesselBlockStride);
            // import connected vessel data
            H5Dread(connectedVesselDatasetId, H5T_NATIVE_LLONG,
                connectedVesselMemSpaceId, connectedVesselDataSpaceId, ftplId,
                connectedVesselIndexArray);
            H5Sclose(connectedVesselMemSpaceId);
            // add connected vessel data to vessels
            long long* startPointer = connectedVesselIndexArray;
            for (long long i = 0; i < maxNumVesselsPerWrite; i++)
            {
                long long vesselIndex = i + remainingNumberOfVesselsToExport + k * maxNumVesselsPerWrite;

                //cout << "hi" << endl;
                long long* endPointer = startPointer + partInfo.numConnectedVesselsArray[vesselIndex];

                vector<long long> temp(startPointer, endPointer);
                std::sort(temp.begin(), temp.end());
                //timsort(temp.begin(), temp.end());

                vessels[vesselIndex].setConnectedVessels(temp);
                startPointer = endPointer;
            }
            connectedVesselStartOffset += connectedVesselCount[0];
        }
        else
        {
            nodeCount[0] = 0;
            connectedVesselCount[0] = 0;
            hid_t nullSpaceId = H5Screate_simple(2, nodeCount, NULL);
            H5Dread(nodeDatasetId, H5T_NATIVE_LLONG, nullSpaceId, nullSpaceId, ftplId, nodeArray);
            H5Sclose(nullSpaceId);
            nullSpaceId = H5Screate_simple(2, connectedVesselCount, NULL);
            H5Dread(connectedVesselDatasetId, H5T_NATIVE_LLONG, nullSpaceId, nullSpaceId, ftplId, connectedVesselIndexArray);
            H5Sclose(nullSpaceId);
        }
    }
    H5Dclose(nodeDatasetId);
    H5Dclose(connectedVesselDatasetId);

    delete[] nodeArray;
    free(connectedVesselIndexArray);
}

void FindAndLabelBcNodes(hid_t fileId, std::vector<EdgeFinderVessel>& vessels)
{
    //double time = -omp_get_wtime();
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t datasetId = OpenHdfDataset(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET);
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t numBcs = 0;
    H5Sget_simple_extent_dims(dataspaceId, &numBcs, NULL);
    long long* bcNodes = new long long[numBcs];

    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    H5Dread(datasetId, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, dtplId, bcNodes);

    H5Sclose(dataspaceId);
    H5Pclose(dtplId);
    H5Dclose(datasetId);
    //if (!mpiRank) cout << "Import BC node time: " << time +omp_get_wtime() << endl;
    const long long numLoops = vessels.size();
    //time = -omp_get_wtime();
//#pragma omp parallel for schedule(guided, 256)
    for (long long i = 0; i < numLoops; ++i)
    {
        EdgeFinderVessel& vessel = vessels[i];
        for (int k = 0; k < numBcs; ++k)
        {
            if (vessel.getStartNode() == bcNodes[k])
            {
                vessels[i].setStartNodeIsBcNode(true);
            }

            if (vessel.getEndNode() == bcNodes[k])
            {
                vessels[i].setEndNodeisBcNode(true);
            }
        }
    }

    delete[] bcNodes;
    //if (!mpiRank) cout << "Label BC nodes time: " << time + omp_get_wtime() << endl;
}

void FindAndLabelBcNodesBcast(hid_t fileId, std::vector<EdgeFinderVessel>& vessels)
{
    //double time = -omp_get_wtime();
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    vector<long long> bcNodes = ImportAndBcastSmallDataset<long long>(fileId, FLOW_GROUP_NAME, BOUNDARY_CONDITION_NODE_DATASET, MPI_LONG_LONG_INT);

    // if (!mpiRank)
    //     cout << "Import BC node time: " << time + omp_get_wtime() << endl;
    const long long numLoops = vessels.size();
//     time = -omp_get_wtime();
// #pragma omp parallel for schedule(guided, 256)
    for (long long i = 0; i < numLoops; ++i)
    {
        EdgeFinderVessel& vessel = vessels[i];
        for (int k = 0; k < bcNodes.size(); ++k)
        {
            if (vessel.getStartNode() == bcNodes[k])
            {
                vessels[i].setStartNodeIsBcNode(true);
            }

            if (vessel.getEndNode() == bcNodes[k])
            {
                vessels[i].setEndNodeisBcNode(true);
            }
        }
    }

    // if (!mpiRank)
    //     cout << "Label BC nodes time: " << time + omp_get_wtime() << endl;
}

void ImportVesselsEdgeFinder(hid_t fileId, PartitionInformation& partInfo, std::vector<EdgeFinderVessel>& vessels)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    ImportNumVesselsInPartition(fileId, mpiRank, partInfo);

    ImportConnectedVesselIndexes(fileId, mpiRank, partInfo);

    vessels.reserve(partInfo.numVesselsInPartition);
   
    ImportDataAndCreateVessels(partInfo, fileId, mpiRank, vessels);
    
    FindAndLabelBcNodes(fileId, vessels);
}

void CreateAndExportEdgeNodeIndexes(hid_t fileId, int totalCount, int startIndex, short localCount, long long* data)
{
    hid_t ftplId = CreateCollectiveDataTransferPropertiesList();
    //create dataset
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;

    //! swapping dims to C style (was on Fortran style?)
    // dims[0] = 1;
    // dims[1] = totalCount;
    dims[0] = totalCount;
    dims[1] = 1;
    
    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += EDGE_NODE_INDEX_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_STD_I64LE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    H5Sclose(dataspaceId);

    //export data
    hid_t memSpaceId;
    hsize_t start[2], count[2], blockStride[2];
    // dims[0] = 1;
    // dims[1] = localCount;
    dims[0] = localCount;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    // count[0] = 1;
    // count[1] = localCount;
    count[0] = localCount;
    count[1] = 1;

    memSpaceId = H5Screate_simple(1, dims, NULL);
    dataspaceId = H5Dget_space(datasetId);

    // start[0] = 0;
    // start[1] = startIndex;
    start[0] = startIndex;
    start[1] = 0;

    // std::cout << std::endl << "in export function" << std::endl;
    // std::cout << data[0] << ", " << data[1] << std::endl;
    // std::cout << "localCount = " << localCount << std::endl;

    //H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, NULL, count, NULL);

    H5Dwrite(datasetId, H5T_NATIVE_LLONG, memSpaceId, dataspaceId, ftplId, data);

    //close stuff
    H5Dclose(datasetId);
    H5Sclose(dataspaceId);
    H5Pclose(ftplId);
}

void CreateAndExportEdgeNodeHostPartition(hid_t fileId, int totalCount, int startIndex, short localCount, short* data)
{
    hid_t ftplId = CreateCollectiveDataTransferPropertiesList();
    //create dataset
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;

    // dims[0] = 1;
    // dims[1] = totalCount;
    dims[0] = totalCount;
    dims[1] = 1;

    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += EDGE_NODE_HOST_PARTITION_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_STD_I16LE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    H5Sclose(dataspaceId);

    //export data
    hid_t memSpaceId;
    hsize_t start[2], count[2], blockStride[2];
    // dims[0] = 1;
    // dims[1] = localCount;
    dims[0] = localCount;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    // count[0] = 1;
    // count[1] = localCount;
    count[0] = localCount;
    count[1] = 1;

    memSpaceId = H5Screate_simple(1, dims, NULL);
    dataspaceId = H5Dget_space(datasetId);

    // start[0] = 0;
    // start[1] = startIndex;
    start[0] = startIndex;
    start[1] = 0;

    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(datasetId, H5T_NATIVE_SHORT, memSpaceId, dataspaceId, ftplId, data);

    //close stuff
    H5Dclose(datasetId);
    H5Sclose(dataspaceId);
    H5Pclose(ftplId);
}

void CreateAndExportEdgeNodeNumAssociatedPartitions(hid_t fileId, int totalCount, int startIndex, short localCount, short* data)
{
    hid_t ftplId = CreateCollectiveDataTransferPropertiesList();
    //create dataset
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;

    // dims[0] = 1;
    // dims[1] = totalCount;
    dims[0] = totalCount;
    dims[1] = 1;

    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += EDGE_NODE_NUM_ASSOCIATED_PARTITIONS_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_STD_I16LE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    H5Sclose(dataspaceId);

    //export data
    hid_t memSpaceId;
    hsize_t start[2], count[2], blockStride[2];
    // dims[0] = 1;
    // dims[1] = localCount;
    dims[0] = localCount;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    // count[0] = 1;
    // count[1] = localCount;
    count[0] = localCount;
    count[1] = 1;

    memSpaceId = H5Screate_simple(1, dims, NULL);
    dataspaceId = H5Dget_space(datasetId);

    // start[0] = 0;
    // start[1] = startIndex;
    start[1] = 0;
    start[0] = startIndex;

    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(datasetId, H5T_NATIVE_SHORT, memSpaceId, dataspaceId, ftplId, data);

    //close stuff
    H5Dclose(datasetId);
    H5Sclose(dataspaceId);
    H5Pclose(ftplId);
}

void CreateAndExportEdgeNodeAssociatedPartitions(hid_t fileId, int totalCount, int startIndex, short localCount, short* data)
{
    hid_t ftplId = CreateCollectiveDataTransferPropertiesList();
    //create dataset
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    std::string datasetPath;

    // dims[0] = 1;
    // dims[1] = totalCount;
    dims[1] = 1;
    dims[0] = totalCount;

    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += EDGE_NODE_ASSOCIATED_PARTITIONS_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_STD_I16LE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    H5Sclose(dataspaceId);

    //export data
    hid_t memSpaceId;
    hsize_t start[2], count[2], blockStride[2];
    // dims[0] = 1;
    // dims[1] = localCount;
    dims[1] = 1;
    dims[0] = localCount;

    blockStride[0] = 1;
    blockStride[1] = 1;
    // count[0] = 1;
    // count[1] = localCount;
    count[1] = 1;
    count[0] = localCount;

    memSpaceId = H5Screate_simple(1, dims, NULL);
    dataspaceId = H5Dget_space(datasetId);

    // start[0] = 0;
    // start[1] = startIndex;
    start[1] = 0;
    start[0] = startIndex;

    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(datasetId, H5T_NATIVE_SHORT, memSpaceId, dataspaceId, ftplId, data);

    //close stuff
    H5Dclose(datasetId);
    H5Sclose(dataspaceId);
    H5Pclose(ftplId);
}

//!!! NEEDS FIXING
void ExportEdgeNodesToFile(hid_t fileId, std::set<EdgeNode>& localEdgeNodes, PartitionInformation& partInfo)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    if (!mpiRank) cout << "exporting Edge Node Data" << endl;

    int totalNumEdgeNodes = 0, numLocalEdgeNodesToExport = 0, totalNumAssociatedPartitions = 0,
        numLocalAssociatedPartitionsToExport = 0, localAssociatedEdgeNodeStart = 0, localEdgeNodeStart = 0;

    for (set<EdgeNode>::iterator itr = localEdgeNodes.begin();
         itr != localEdgeNodes.end(); itr++)
    {
        if (itr->getHostPartition() == partInfo.partitionNum)
        {
            numLocalEdgeNodesToExport++;
            numLocalAssociatedPartitionsToExport += ((itr->getConnectedPartitions().size()) - 1);
        }
    }

    MPI_Allreduce(&numLocalEdgeNodesToExport, &totalNumEdgeNodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Exscan(&numLocalEdgeNodesToExport, &localEdgeNodeStart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (!mpiRank) {localEdgeNodeStart = 0;}

    MPI_Allreduce(&numLocalAssociatedPartitionsToExport, &totalNumAssociatedPartitions, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Exscan(&numLocalAssociatedPartitionsToExport, &localAssociatedEdgeNodeStart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (!mpiRank) {localAssociatedEdgeNodeStart = 0;}

    //prepare data
    long long* index = new long long[numLocalEdgeNodesToExport];
    short* hostPartition = new short[numLocalEdgeNodesToExport];
    short* numAssociatedParts = new short[numLocalEdgeNodesToExport];
    short* associatedParts = new short[numLocalAssociatedPartitionsToExport];

    short nodeIndex = 0, associatedPartIndex = 0;
    for (set<EdgeNode>::iterator itr = localEdgeNodes.begin(); itr != localEdgeNodes.end(); itr++)
    {
        if (itr->getHostPartition() == partInfo.partitionNum)
        {
            index[nodeIndex] = itr->getNode();
            //std::cout << "partInfo.partitionNum = " << partInfo.partitionNum << std::endl;
            //std::cout << "Edge Node Index: " << index[nodeIndex] << ", on rank " << mpiRank << std::endl;

            hostPartition[nodeIndex] = itr->getHostPartition();
            //std::cout << "Host Partition: " << hostPartition[nodeIndex] << ", on rank " << mpiRank << std::endl;

            auto currentAssociatedPartitionSet = itr->getConnectedPartitions();
            numAssociatedParts[nodeIndex] = (currentAssociatedPartitionSet.size() - 1);
            for (auto& partition : currentAssociatedPartitionSet)
            {
                if (partition != itr->getHostPartition())
                {
                    associatedParts[associatedPartIndex++] = partition;
                }
            }
            nodeIndex++;
        }
    }

    // std::cout << index[0] << ", " << index[1] << ", on rank " << mpiRank << std::endl;
    // std::cout << "totalNumEdgeNodes = " << totalNumEdgeNodes << ", on rank " << mpiRank << std::endl;
    // std::cout << "localEdgeNodeStart = " << localEdgeNodeStart << ", on rank " << mpiRank << std::endl;
    // std::cout << "numLocalEdgeNodesToExport = " << numLocalEdgeNodesToExport << ", on rank " << mpiRank << std::endl;

    //export node Index
    //! double check - only showing one edge node
    CreateAndExportEdgeNodeIndexes(fileId, totalNumEdgeNodes, localEdgeNodeStart, numLocalEdgeNodesToExport, index);

    //export host partition
    //! double check
    CreateAndExportEdgeNodeHostPartition(fileId, totalNumEdgeNodes, localEdgeNodeStart, numLocalEdgeNodesToExport, hostPartition);

    //export num associated parts
    //! double check
    CreateAndExportEdgeNodeNumAssociatedPartitions(fileId, totalNumEdgeNodes, localEdgeNodeStart, numLocalEdgeNodesToExport, numAssociatedParts);

    //export associated parts
    //! double check
    CreateAndExportEdgeNodeAssociatedPartitions(fileId,
        totalNumAssociatedPartitions, localAssociatedEdgeNodeStart,
        numLocalAssociatedPartitionsToExport, associatedParts);

    delete[] index;
    delete[] hostPartition;
    delete[] numAssociatedParts;
    delete[] associatedParts;
}