#include <iostream>
#include <new>
#include <numeric>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

//#include <vascular/core/Geometry.h>
//#include <vascular/core/GlobalDefs.h>
#include "VesselGenerator/coordinate.h"
//#include "VesselGenerator/vessel.h"
#include "VesselGenerator/ImportExportCommon.h"
#include "PreprocessorVessel.h"
#include "SortedVesselExporter.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

//typedef PreprocessorVessel Vessel;
typedef vector<PreprocessorVessel> VesselVector;
typedef unsigned long ulong;
//typedef Point3D<double> Point;

#define IDX(i, j, n) n* i + j

void CreateGeomExportVector(long long numVesselToProcess,
    long long localVesselOffset, VesselVector& vessels,
    vector<double>& dataArray)
{
    dataArray.resize(numVesselToProcess * 7);
    for (long long j = 0; j < numVesselToProcess; j++)
    {
        long long vesselIndex = localVesselOffset + j;
        size_t dataTableIndex = j * 7;
        Coordinate start = vessels[vesselIndex].get_startingPoint();
        Coordinate end = vessels[vesselIndex].get_endingPoint();
        dataArray[dataTableIndex++] = start.x;
        dataArray[dataTableIndex++] = start.y;
        dataArray[dataTableIndex++] = start.z;
        dataArray[dataTableIndex++] = end.x;
        dataArray[dataTableIndex++] = end.y;
        dataArray[dataTableIndex++] = end.z;
        dataArray[dataTableIndex++] = vessels[vesselIndex].get_radius();
    }
}

void PerfromBlockGeomExport(hid_t datasetId, long long numVesselsToExport,
    long long fileVesselOffset, vector<double>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    
    if (remainingNumberOfVesselsToExport > 0)
        CollectiveHdfBlockExport<double, 2, 7>(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray);
    else
        CollectiveHdfNullExport<double, 2, 7>(datasetId);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) * 7;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            CollectiveHdfBlockExport<double, 2, 7>(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
        else
        {
            CollectiveHdfNullExport<double, 2, 7>(datasetId);
        }
    }
}

void PerfromBlockGeomExportIndependent(hid_t datasetId,
    long long numVesselsToExport, long long fileVesselOffset,
    vector<double>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    
    if (remainingNumberOfVesselsToExport > 0)
        IndependentHdfBlockExport<double, 2, 7>(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) * 7;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            IndependentHdfBlockExport<double, 2, 7>(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
    }
}

void ExportGeomData(hid_t fileId, VesselVector vessels, hsize_t* vesselOffsets,
    long long* localNumVesselsPerPartition, int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    vector<double> dataArray;
    hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, GEOM_DATASET);
    for (int i = 0; i < numParts; i++)
    {
        const int partitionToExport = (i + mpiRank) % mpiSize;
        const long long localVesselOffset = std::accumulate(localNumVesselsPerPartition, localNumVesselsPerPartition + partitionToExport, 0LL);
        CreateGeomExportVector(localNumVesselsPerPartition[partitionToExport], localVesselOffset, vessels, dataArray);
        PerfromBlockGeomExport(datasetId, localNumVesselsPerPartition[partitionToExport], vesselOffsets[partitionToExport], dataArray);
    }
    H5Dclose(datasetId);
}

hid_t CreatOriginalIndexDataset(hid_t fileId, hsize_t numVessels, std::string groupName, std::string datasetName)
{
    hid_t dataspaceId, datasetId = -1;
    herr_t status;
    hid_t dcpl;
    hsize_t dims[2];
    string datasetPath = groupName + datasetName;

    H5E_auto1_t origFunction = nullptr;
    void** client_data = nullptr;
    H5Eget_auto1(&origFunction, client_data);
    H5Eset_auto1(NULL, NULL);
    datasetId = H5Dopen2(fileId, datasetPath.c_str(), H5P_DEFAULT);
    H5Eset_auto1(origFunction, client_data);
    if (datasetId > 0)
    {
        return datasetId;
    }

    dims[0] = numVessels;
    dims[1] = 1;
    hsize_t maxDims[2] = {H5S_UNLIMITED};
    dataspaceId = H5Screate_simple(1, dims, NULL);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[2];
    if (dims[0] > 524286)
    {
        chunk[0] = dims[0] / 1024;
        chunk[1] = 1;
        //status = H5Pset_chunk(dcpl, 1, chunk);
    }
    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_NATIVE_LLONG,
        dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(dataspaceId);

    return datasetId;
}

void ExportOriginalIndexes(hid_t fileId, VesselVector vessels,
    hsize_t* vesselOffsets, long long* localNumVesselsPerPartition, int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long numLocalVessels = vessels.size();
    long long totalNumVessels = 0;
    MPI_Allreduce(&numLocalVessels, &totalNumVessels, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    hid_t datasetId = CreatOriginalIndexDataset(fileId, totalNumVessels, GEOM_GROUP_NAME, ORGINAL_INDEX_DATASET);

    vector<long long> dataArray;
    hsize_t memDim = vessels.size();
    hsize_t localVesselOffset = 0;
    for (int i = 0; i < numParts; i++) {    //set dataspace
        const int partitionToExport = (i + mpiRank) % mpiSize;
        dataArray.resize(localNumVesselsPerPartition[partitionToExport]);
        localVesselOffset = std::accumulate(localNumVesselsPerPartition, localNumVesselsPerPartition + partitionToExport, 0LL);
        
        for (long long j = 0; j < localNumVesselsPerPartition[partitionToExport]; ++j) {
            dataArray[j] = vessels[j + localVesselOffset].getInitialArrayPosition();
        }

        CollectiveHdfBlockExport<long long>(datasetId, vesselOffsets[partitionToExport],
            localNumVesselsPerPartition[partitionToExport], dataArray);
    }
    H5Dclose(datasetId);
}

void ExportOriginalIndexesOld(hid_t fileId, VesselVector vessels,
    hsize_t* vesselOffsets, long long* localNumVesselsPerPartition, int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long numLocalVessels = vessels.size();
    long long totalNumVessels = 0;
    MPI_Allreduce(&numLocalVessels, &totalNumVessels, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    //cout << mpiRank << ": Entering Dataset Creator " << endl;
    hid_t datasetId = CreatOriginalIndexDataset(fileId, totalNumVessels, GEOM_GROUP_NAME, ORGINAL_INDEX_DATASET);

    //cout << mpiRank << ": Exiting Dataset Creator " << endl;
    unsigned long long* dataArray = new unsigned long long[vessels.size()];

    size_t dataTableIndex = 0;
    for (VesselVector::const_iterator itr = vessels.begin(); itr != vessels.end(); itr++) {
        dataArray[dataTableIndex++] = itr->getInitialArrayPosition();
    }

    hid_t dataTransferPropertiesList = H5Pcreate(H5P_DATASET_XFER);
    if (mpiSize > 1)
        H5Pset_dxpl_mpio(dataTransferPropertiesList, H5FD_MPIO_COLLECTIVE);
    hsize_t memDim = vessels.size();
    hsize_t localMemOffset = 0;
    for (int i = 0; i < numParts; i++) {    //set dataspace
        hid_t dataspaceId = H5Dget_space(datasetId);

        //set initial hyperslab
        hsize_t start[1], block[1], count[1], stride[1];

        stride[0] = 1;
        //stride[1] = 1;
        block[0] = 1;
        //block[1] = 1;
        start[0] = vesselOffsets[i];
        //start[1] = 0;
        count[0] = localNumVesselsPerPartition[i];
        //count[1] = 1;
        //cout << "Partition: " << i << " MPI Rank: " << mpiRank << " Start: " << start[0] << " Count: " << count[0] << endl;
        H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, stride, count, block);

        //use HDF5 write
        hid_t memSpaceId = H5Screate_simple(1, &memDim, NULL);
        H5Sselect_hyperslab(memSpaceId, H5S_SELECT_SET, &localMemOffset, stride, count, block);
        H5Dwrite(datasetId, H5T_NATIVE_ULLONG, memSpaceId, dataspaceId, dataTransferPropertiesList, dataArray);

        H5Sclose(dataspaceId);
        H5Sclose(memSpaceId);
        localMemOffset += localNumVesselsPerPartition[i];
    }
    H5Dclose(datasetId);
    H5Pclose(dataTransferPropertiesList);
    delete[] dataArray;
}

void PerfromBlockNodeExport(hid_t datasetId, long long numVesselsToExport, long long fileVesselOffset, vector<long long>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        CollectiveHdfBlockExport<long long, 2, 2>(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray);
    else
        CollectiveHdfNullExport<long long, 2, 2>(datasetId);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) * 2;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            CollectiveHdfBlockExport<long long, 2, 2>(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
        else
        {
            CollectiveHdfNullExport<long long, 2, 2>(datasetId);
        }
    }
}

void PerfromBlockNodeExportIndependent(hid_t datasetId, long long numVesselsToExport, long long fileVesselOffset, vector<long long>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        IndependentHdfBlockExport<long long, 2, 2>(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) * 2;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            IndependentHdfBlockExport<long long, 2, 2>(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
    }
}

void CreateNodeExportVector(long long numVesselToProcess,
    long long localVesselOffset, VesselVector& vessels,
    vector<long long>& dataArray)
{
    dataArray.resize(numVesselToProcess * 7);
    for (long long j = 0; j < numVesselToProcess; j++)
    {
        long long vesselIndex = localVesselOffset + j;
        size_t dataTableIndex = j * 2;
        dataArray[dataTableIndex++] = vessels[vesselIndex].getNode1();
        dataArray[dataTableIndex++] = vessels[vesselIndex].getNode2();
    }
}

void ExportNodeData(hid_t fileId, VesselVector vessels, hsize_t* vesselOffsets,
    long long* localNumVesselsPerPartition, int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    vector<long long> dataArray;
    hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);
    for (int i = 0; i < numParts; i++)
    {
        const int partitionToExport = (i + mpiRank) % mpiSize;
        long long localVesselOffset = std::accumulate(localNumVesselsPerPartition,
                localNumVesselsPerPartition + partitionToExport, 0LL);
        CreateNodeExportVector(localNumVesselsPerPartition[partitionToExport], localVesselOffset, vessels, dataArray);

        PerfromBlockNodeExport(datasetId, localNumVesselsPerPartition[partitionToExport],
            vesselOffsets[partitionToExport], dataArray);
        //PerfromBlockNodeExportIndependent(datasetId, localNumVesselsPerPartition[i], vesselOffsets[i], dataArray);
    }
    H5Dclose(datasetId);
}

void ExportPartitionAttributes(hid_t fileId, int numParts)
{
    //open group
    string attributeName = "Num_Partitions";
    hsize_t dims = 1;
    hid_t dataSpaceId = H5Screate_simple(1, &dims, NULL);

    hid_t attributeId = H5Aopen(fileId, attributeName.c_str(), H5P_DEFAULT);

    H5Awrite(attributeId, H5T_NATIVE_INT, &numParts);

    H5Sclose(dataSpaceId);
    H5Aclose(attributeId);

    attributeName = "Is_Partitioned";
    dataSpaceId = H5Screate_simple(1, &dims, NULL);

    int wasPartitioned = 1;

    attributeId = H5Aopen(fileId, attributeName.c_str(), H5P_DEFAULT);

    H5Awrite(attributeId, H5T_NATIVE_INT, &wasPartitioned);

    H5Sclose(dataSpaceId);
    H5Aclose(attributeId);
}

void PerformBlockNumConnectedVesselExport(hid_t datasetId,
    long long numVesselsToExport, long long fileVesselOffset,
    vector<short>& dataArray, long long partitionDataOffset)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        CollectiveHdfBlockExport(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray, partitionDataOffset);
    else
        CollectiveHdfNullExport<short>(datasetId);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) + partitionDataOffset;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            CollectiveHdfBlockExport(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
        else
        {
            CollectiveHdfNullExport<short>(datasetId);
        }
    }
}

void PerformBlockNumConnectedVesselExportIndependent(hid_t datasetId,
    long long numVesselsToExport, long long fileVesselOffset,
    vector<short>& dataArray, long long partitionDataOffset)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = numVesselsToExport % maxNumVesselsPerWrite;
    int numLoopsNeeded = numVesselsToExport / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);
    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        IndependentHdfBlockExport(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray, partitionDataOffset);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite) + partitionDataOffset;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
            IndependentHdfBlockExport(datasetId, fileOffset, maxNumVesselsPerWrite, dataArray, dataOffset);
        }
    }
}

void PerformBlockConnectedVesselExport(hid_t datasetId,
    long long totalNumElementsToExport, long long fileVesselOffset,
    long long partitionDataOffset, vector<long long>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    long long maxNumberOfConnectedVesselsToWrite = maxNumVesselsPerWrite * 4;
    long long remainingNumberOfVesselsToExport = totalNumElementsToExport % maxNumberOfConnectedVesselsToWrite;
    int numLoopsNeeded = totalNumElementsToExport / maxNumberOfConnectedVesselsToWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        CollectiveHdfBlockExport(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray, partitionDataOffset);
    else
        CollectiveHdfNullExport<long long>(datasetId);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumberOfConnectedVesselsToWrite) + partitionDataOffset;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumberOfConnectedVesselsToWrite;
            CollectiveHdfBlockExport(datasetId, fileOffset, maxNumberOfConnectedVesselsToWrite, dataArray, dataOffset);
        }
        else
        {
            CollectiveHdfNullExport<long long>(datasetId);
        }
    }
}

void PerformBlockConnectedVesselExportIndependent(hid_t datasetId,
    long long totalNumElementsToExport, long long fileVesselOffset,
    long long partitionDataOffset, vector<long long>& dataArray)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    long long maxNumberOfConnectedVesselsToWrite = maxNumVesselsPerWrite * 4;
    long long remainingNumberOfVesselsToExport = totalNumElementsToExport % maxNumberOfConnectedVesselsToWrite;
    int numLoopsNeeded = totalNumElementsToExport / maxNumberOfConnectedVesselsToWrite;
    int globalNumLoopsNeeded = 0;
    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    //cout << "remainder: " << remainingNumberOfVesselsToExport << " Loops: " << numLoopsNeeded << " global loops " << endl;
    if (remainingNumberOfVesselsToExport > 0)
        IndependentHdfBlockExport(datasetId, fileVesselOffset, remainingNumberOfVesselsToExport, dataArray, partitionDataOffset);

    //perform loops
    for (int i = 0; i < globalNumLoopsNeeded; ++i)
    {
        if (i < numLoopsNeeded)
        {
            long long dataOffset = (remainingNumberOfVesselsToExport + i * maxNumberOfConnectedVesselsToWrite) + partitionDataOffset;
            long long fileOffset = fileVesselOffset + remainingNumberOfVesselsToExport + i * maxNumberOfConnectedVesselsToWrite;
            IndependentHdfBlockExport(datasetId, fileOffset, maxNumberOfConnectedVesselsToWrite, dataArray, dataOffset);
        }
    }
}

void ExportConnectedVessels(hid_t fileId, VesselVector vessels,
    hsize_t* vesselOffsets, long long* localNumVesselsPerPartition,
    int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t numVesselsDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NUM_CONNECTED_VESSELS_DATASET);
    hid_t connectedVesselsDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);
    long long numConnectedVessels = 0;

    vector<short> numVesselsDataArray(vessels.size());
    vector<hsize_t> localConnectedVesselCount(numParts);

    size_t numVesselsdataTableIndex = 0;

    for (VesselVector::const_iterator itr = vessels.begin();
         itr != vessels.end(); itr++)
    {
        int vesselcount = itr->getNumConnectedVessels();
        numVesselsDataArray[numVesselsdataTableIndex++] = vesselcount;
        numConnectedVessels += vesselcount;
        localConnectedVesselCount[itr->getPartitionNumber()] += vesselcount;
    }
    vector<hsize_t> vesselDataspaceOffsets(numParts);
    vector<long long> connectedVesselArray(numConnectedVessels);

    //calculate memory offset
    vector<hsize_t> memoryOffset(numParts);
    memoryOffset[0] = 0;
    for (int i = 1; i < numParts; i++)
        memoryOffset[i] = localConnectedVesselCount[i - 1] + memoryOffset[i - 1];

    CalculateConnectedVesselOffsets(vessels, numParts, vesselDataspaceOffsets);

    unsigned long long tableIndex = 0;
    VesselVector::const_iterator arteryItr;
    for (arteryItr = vessels.begin(); arteryItr != vessels.end(); arteryItr++)
    {
        const vector<long long>& connectedVessels = arteryItr->getConnectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin(); itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }

    for (int i = 0; i < numParts; i++)
    {    //set dataspace

        const int partitionToExport = (i + mpiRank) % mpiSize;
        long long partitionDataOffset =
            std::accumulate(localNumVesselsPerPartition,
                localNumVesselsPerPartition + partitionToExport, 0LL);
        //PerformBlockNumConnectedVesselExportIndependent(numVesselsDatasetId,localNumVesselsPerPartition[i],vesselOffsets[i],numVesselsDataArray,partitionDataOffset);

        //PerformBlockConnectedVesselExportIndependent(connectedVesselsDatasetId,localConnectedVesselCount[i],vesselDataspaceOffsets[i],memoryOffset[i],connectedVesselArray);

        PerformBlockNumConnectedVesselExport(numVesselsDatasetId,
            localNumVesselsPerPartition[partitionToExport],
            vesselOffsets[partitionToExport], numVesselsDataArray,
            partitionDataOffset);

        PerformBlockConnectedVesselExport(connectedVesselsDatasetId,
            localConnectedVesselCount[partitionToExport],
            vesselDataspaceOffsets[partitionToExport],
            memoryOffset[partitionToExport], connectedVesselArray);
    }

    H5Dclose(numVesselsDatasetId);
    H5Dclose(connectedVesselsDatasetId);
}

void ExportVesselsPerPartition(hid_t fileId, long long* localNumVesselsPerPartition, int numParts)
{
    string groupName = "/";

    hid_t dataspaceId, datasetId = -1;
    herr_t status;
    hid_t dcpl;
    hsize_t dims[2];
    string datasetPath = groupName + NUM_VESSELS_PER_PARTITION_DATASET;

    H5E_auto1_t origFunction = nullptr;
    void** client_data = nullptr;
    H5Eget_auto1(&origFunction, client_data);
    H5Eset_auto1(NULL, NULL);
    datasetId = H5Dopen2(fileId, datasetPath.c_str(), H5P_DEFAULT);
    H5Eset_auto1(origFunction, client_data);
    if (datasetId < 0)
    {
        dims[0] = numParts;
        dims[1] = 1;
        hsize_t maxDims[2] = {H5S_UNLIMITED};
        dataspaceId = H5Screate_simple(1, dims, NULL);

        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk[2];
        if (dims[0] > 32)
        {
            chunk[0] = dims[0] / 16;
            chunk[1] = 1;
            status = H5Pset_chunk(dcpl, 1, chunk);
        }
        datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_NATIVE_LLONG,
            dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        status = H5Pclose(dcpl);
        status = H5Sclose(dataspaceId);
    }

    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long* numVesselsPerPartition = new long long[numParts]();

    MPI_Reduce(localNumVesselsPerPartition, numVesselsPerPartition, numParts,
        MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpiRank == 0)
        H5Dwrite(datasetId, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            numVesselsPerPartition);

    delete[] numVesselsPerPartition;
    H5Dclose(datasetId);
}

void ExportVesselDataSorted(hid_t fileId, std::vector<PreprocessorVessel>& vessels, int numParts)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    //Count number of vessels in each partition
    long long* localNumVesselsPerPartition = new long long[numParts]();

    for (VesselVector::const_iterator itr = vessels.begin();
         itr != vessels.end(); itr++)
    {
        localNumVesselsPerPartition[itr->getPartitionNumber()]++;
    }

    ExportPartitionAttributes(fileId, numParts);
    hsize_t* vesselOffsets = new hsize_t[numParts]();

    CalculateVesselOffsets(vessels, numParts, localNumVesselsPerPartition, vesselOffsets);

    ExportGeomData(fileId, vessels, vesselOffsets, localNumVesselsPerPartition, numParts);

    //Export Node Data
    ExportNodeData(fileId, vessels, vesselOffsets, localNumVesselsPerPartition, numParts);
    
    ExportConnectedVessels(fileId, vessels, vesselOffsets, localNumVesselsPerPartition, numParts);
    
    ExportVesselsPerPartition(fileId, localNumVesselsPerPartition, numParts);
    
    H5Fflush(fileId, H5F_SCOPE_LOCAL);
    
    delete[] localNumVesselsPerPartition;
    delete[] vesselOffsets;
}



void CalculateVesselOffsets(std::vector<PreprocessorVessel>& vessels,
    int numParts, long long*& localNumVesselsPerPartition,
    hsize_t*& outputOffsets)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    //Communicate info to all nodes
    long long* numVesselsPerProcPerPartition = new long long[numParts * mpiSize]();

    MPI_Allgather(localNumVesselsPerPartition, numParts, MPI_LONG_LONG_INT,
        numVesselsPerProcPerPartition, numParts, MPI_LONG_LONG_INT, MPI_COMM_WORLD);


    hsize_t* partitionOffset = new hsize_t[numParts]();
    partitionOffset[0] = 0;
    for (int i = 1; i < numParts; i++)
    {
        partitionOffset[i] = partitionOffset[i - 1];
        for (int j = 0; j < mpiSize; j++)
        {
            partitionOffset[i] += numVesselsPerProcPerPartition[IDX(j, i - 1, numParts)];
        }
    }

    if (outputOffsets == nullptr)
        outputOffsets = new hsize_t[numParts]();

    for (int i = 0; i < numParts; i++)
    {
        outputOffsets[i] = partitionOffset[i];
        for (int j = 0; j < mpiRank; j++)
            outputOffsets[i] += numVesselsPerProcPerPartition[IDX(j, i, numParts)];
    }

    delete[] numVesselsPerProcPerPartition;
    delete[] partitionOffset;
}

// void CalculateConnectedVesselOffsets(std::vector<PreprocessorVessel>& vessels, int numParts, hsize_t*& outputOffsets)
// {
//     int mpiSize, mpiRank;

//     MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

//     long long* localNumVesselsPerPartition = new long long[numParts]();

//     for (VesselVector::const_iterator itr = vessels.begin();
//          itr != vessels.end(); itr++)
//     {
//         localNumVesselsPerPartition[itr->getPartitionNumber()] +=
//             itr->getNumConnectedVessels();
//     }

//     //cout << mpiRank << ": Performing Communication " << endl;
//     long long* numVesselsPerProcPerPartition =
//         new long long[numParts * mpiSize]();

//     MPI_Allgather(localNumVesselsPerPartition, numParts, MPI_LONG_LONG,
//         numVesselsPerProcPerPartition, numParts, MPI_LONG_LONG, MPI_COMM_WORLD);
//     //cout << mpiRank << ": calculating partition offset " << endl;
//     //outputOffsets = new hsize_t[numParts];
//     hsize_t* partitionOffset = new hsize_t[numParts];
//     partitionOffset[0] = 0;

//     for (int i = 1; i < numParts; i++)
//     {
//         partitionOffset[i] = partitionOffset[i - 1];
//         for (int j = 0; j < mpiSize; j++)
//             partitionOffset[i] +=
//                 numVesselsPerProcPerPartition[IDX(j, i - 1, numParts)];
//     }

//     if (outputOffsets == nullptr)
//         outputOffsets = new hsize_t[numParts]();

//     for (int i = 0; i < numParts; i++)
//     {
//         outputOffsets[i] = partitionOffset[i];
//         for (int j = 0; j < mpiRank; j++)
//             outputOffsets[i] +=
//                 numVesselsPerProcPerPartition[IDX(j, i, numParts)];
//     }

//     delete[] localNumVesselsPerPartition;
//     delete[] numVesselsPerProcPerPartition;
//     delete[] partitionOffset;
// }

void CalculateConnectedVesselOffsets(std::vector<PreprocessorVessel>& vessels, int numParts, std::vector<hsize_t>& outputOffsets)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long* localNumVesselsPerPartition = new long long[numParts]();

    for (VesselVector::const_iterator itr = vessels.begin(); itr != vessels.end(); itr++)
    {
        localNumVesselsPerPartition[itr->getPartitionNumber()] += itr->getNumConnectedVessels();
    }

    //cout << mpiRank << ": Performing Communication " << endl;
    long long* numVesselsPerProcPerPartition = new long long[numParts * mpiSize]();

    MPI_Allgather(localNumVesselsPerPartition, numParts, MPI_LONG_LONG,
        numVesselsPerProcPerPartition, numParts, MPI_LONG_LONG, MPI_COMM_WORLD);
    //cout << mpiRank << ": calculating partition offset " << endl;
    //outputOffsets = new hsize_t[numParts];
    hsize_t* partitionOffset = new hsize_t[numParts];
    partitionOffset[0] = 0;

    for (int i = 1; i < numParts; i++)
    {
        partitionOffset[i] = partitionOffset[i - 1];
        for (int j = 0; j < mpiSize; j++)
            partitionOffset[i] += numVesselsPerProcPerPartition[IDX(j, i - 1, numParts)];
    }

    outputOffsets.resize(numParts);
    for (int i = 0; i < numParts; i++)
    {
        outputOffsets[i] = partitionOffset[i];
        for (int j = 0; j < mpiRank; j++)
            outputOffsets[i] += numVesselsPerProcPerPartition[IDX(j, i, numParts)];
    }

    delete[] localNumVesselsPerPartition;
    delete[] numVesselsPerProcPerPartition;
    delete[] partitionOffset;
}