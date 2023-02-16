#include <bitset>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "ImportExportCommon.h"
#include "GeneratorSharedHdfFunctions.h"
#include "vessel.h"
#include "VesselNetworkExporter.h"
#include "coordinate.h"

using std::bitset;
using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef vector<Vessel> VesselVector;
typedef unsigned long ulong;

#define NUM_ITEMS_PER_CHUNK 1024

void WriteSubVesselsToFile(long long subArterialTreeNumberOfVessels,
    long long rootArteriealTreeNumberOfVessels, hid_t file_id,
    VesselVector& arteries,
    VesselVector& veins)
{
    long long totalNumberOfVesselsToExport = 2 * subArterialTreeNumberOfVessels;
    long long totalNumberOfRootVessels = 2 * rootArteriealTreeNumberOfVessels;
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    // cout << "Create Geom Dataset" << endl;
    hid_t geomDatasetId = OpenHdfDataset(file_id, GEOM_GROUP_NAME, GEOM_DATASET);

    vector<double> geomVector(subArterialTreeNumberOfVessels * 7);
    
    ConvertVesselVectorsToGeometryArray(subArterialTreeNumberOfVessels, arteries, geomVector.data());

    // Determine number of blocks for successful write
    const long long remainingNumberOfVesselsToExport =
        subArterialTreeNumberOfVessels % maxNumVesselsPerWrite;
    const int numLoopsNeeded =
        subArterialTreeNumberOfVessels / maxNumVesselsPerWrite;

    // cout << numLoopsNeeded << ", " << remainingNumberOfVesselsToExport<<endl;
    size_t memoryOffset;

    // hid_t dtpl = CreateCollectiveDataTransferPropertiesList();
    // hid_t memSpaceId, dataSpaceId;
    hsize_t dims[2] = {0};
    hsize_t start[2] = {0}, count[2] = {0}, blockStride[2] = {0};
    dims[0] = remainingNumberOfVesselsToExport;
    dims[1] = 7;
    // memSpaceId = H5Screate_simple(2, dims, NULL);
    // dataSpaceId = H5Dget_space(geomDatasetId);

    blockStride[0] = 1;
    blockStride[1] = 1;

    count[0] = remainingNumberOfVesselsToExport;
    count[1] = 7;
    start[1] = 0;

    if (remainingNumberOfVesselsToExport > 0)
    {
        count[0] = remainingNumberOfVesselsToExport;

        start[0] = totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport;

        CollectiveHdfBlockExport<double, 2, 7>(geomDatasetId, start[0], count[0], geomVector);
    }

    count[0] = maxNumVesselsPerWrite;

    for (int i = 0; i < numLoopsNeeded; i++)
    {
        memoryOffset = i * ((count[0]) * count[1]) +
            remainingNumberOfVesselsToExport * count[1];
        start[0] = totalNumberOfRootVessels +
            mpiRank * totalNumberOfVesselsToExport + i * count[0] +
            remainingNumberOfVesselsToExport;

        CollectiveHdfBlockExport<double, 2, 7>(
            geomDatasetId, start[0], count[0], geomVector, memoryOffset);
    }

    // perform write for veins
    ConvertVesselVectorsToGeometryArray(subArterialTreeNumberOfVessels, veins, geomVector.data());
    int numWritesNeeded = 1;

    dims[0] = subArterialTreeNumberOfVessels / numWritesNeeded;
    dims[1] = 7;
    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = totalNumberOfRootVessels +
        mpiRank * totalNumberOfVesselsToExport + subArterialTreeNumberOfVessels;
    start[1] = 0;
    count[0] = subArterialTreeNumberOfVessels / numWritesNeeded;
    count[1] = 7;

    // H5Sselect_hyperslab(dataSpaceId, H5S_SELECT_SET, start, blockStride, count,
    // blockStride);
    // H5Dwrite(geomDatasetId, H5T_NATIVE_DOUBLE, memSpaceId, dataSpaceId, dtpl,
    // geomArray);

    count[0] = remainingNumberOfVesselsToExport;
    count[1] = 7;
    start[1] = 0;

    if (remainingNumberOfVesselsToExport > 0)
    {
        count[0] = remainingNumberOfVesselsToExport;

        start[0] = totalNumberOfRootVessels +
            mpiRank * totalNumberOfVesselsToExport +
            subArterialTreeNumberOfVessels;
        
        CollectiveHdfBlockExport<double, 2, 7>(geomDatasetId, start[0], count[0], geomVector);
    }

    count[0] = maxNumVesselsPerWrite;

    for (int i = 0; i < numLoopsNeeded; i++)
    {
        memoryOffset = i * ((count[0]) * count[1]) +
            remainingNumberOfVesselsToExport * count[1];
        start[0] = totalNumberOfRootVessels +
            mpiRank * totalNumberOfVesselsToExport + i * count[0] +
            remainingNumberOfVesselsToExport + subArterialTreeNumberOfVessels;

        CollectiveHdfBlockExport<double, 2, 7>(geomDatasetId, start[0], count[0], geomVector, memoryOffset);
    }

    H5Dclose(geomDatasetId);
    
    geomVector.clear();

    // add sub nodes to file
    // cout << "Create Node Dataset" << endl;
    hid_t nodeDatasetId = OpenHdfDataset(file_id, GEOM_GROUP_NAME, NODE_DATASET);

    vector<long long> nodeVector(totalNumberOfVesselsToExport * 2);
    ConvertVesselVectorsToNodeArray(totalNumberOfVesselsToExport / 2, arteries, veins, nodeVector.data());

    CollectiveHdfBlockExport<long long, 2, 2>(nodeDatasetId,
        totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport,
        totalNumberOfVesselsToExport, nodeVector);

    H5Dclose(nodeDatasetId);

    // delete[] nodeArray;
    nodeVector.clear();

    //! write sub conductances to matrix
    hid_t conductanceDatasetId = OpenHdfDataset(file_id, FLOW_GROUP_NAME, HEALTHY_CONDUCTANCE_DATASET);

    vector<double> conductanceVector(totalNumberOfVesselsToExport);
    ConvertVesselVectorsToConductanceArray(totalNumberOfVesselsToExport / 2, arteries, veins, conductanceVector.data());
    // for(int i = 0; i < conductanceVector.size(); i++) {
    //     std::cout << conductanceVector[i] << std::endl;
    // }
    CollectiveHdfBlockExport<double, 2, 1>(conductanceDatasetId,
        totalNumberOfRootVessels + mpiRank * totalNumberOfVesselsToExport,
        totalNumberOfVesselsToExport, conductanceVector);

    H5Dclose(conductanceDatasetId);

    conductanceVector.clear();
}

void WriteSubConnectedVesselDataToFile(long numVesselsExport,
    long long rootNumVessels, long long totalNumConnectedVessels,
    hid_t numConnectionsDataset, short* numConnectedVesselsArray,
    hid_t connectionsArrayDataset, long long* connectedVesselArray)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    hid_t dtpl = CreateCollectiveDataTransferPropertiesList();

    hid_t memSpaceId, dataSpaceId;
    hsize_t dims[2];
    hsize_t start[2], count[2], blockStride[2];
    dims[0] = numVesselsExport;
    dims[1] = 1;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(numConnectionsDataset);

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = rootNumVessels + mpiRank * numVesselsExport;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 1;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    // cout <<  << " " << rootNumVessels << endl;
    // cout << "Write Num Connected Vessels" << endl;

    H5Dwrite(numConnectionsDataset, H5T_NATIVE_SHORT, memSpaceId, dataSpaceId,
        dtpl, numConnectedVesselsArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);
    H5Pclose(dtpl);

    H5Dclose(numConnectionsDataset);

    // dtpl = CreateCollectiveDataTransferPropertiesList();
    // const long maxNumVesselsPerWrite = 67108864;
    // cout << subArterialTreeNumberOfVessels << ", " << arteries.size() << ", "
    // << maxNumVesselsPerWrite<< endl;

    const long long remainingNumberOfVesselsToExport =
        totalNumConnectedVessels % maxNumVesselsPerWrite;
    const int numLoopsNeeded = totalNumConnectedVessels / maxNumVesselsPerWrite;
    size_t memoryOffset;

    dims[0] = totalNumConnectedVessels;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank;
    start[1] = 0;
    count[0] = totalNumConnectedVessels;
    count[1] = 1;
    dataSpaceId = H5Dget_space(connectionsArrayDataset);
    if (remainingNumberOfVesselsToExport > 0)
    {
        count[0] = remainingNumberOfVesselsToExport;

        start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank;
        // start[1] = 0;

        memSpaceId = H5Screate_simple(2, count, NULL);
        // cout << "Rank" << mpiRank << ": " << start[0] << ", " << count[0] <<
        // endl;
        H5Sselect_hyperslab(dataSpaceId, H5S_SELECT_SET, start, blockStride,
            count, blockStride);
        // cout << "write" << endl;
        H5Dwrite(connectionsArrayDataset, H5T_NATIVE_LLONG, memSpaceId,
            dataSpaceId, dtpl, connectedVesselArray);
        // cout << "end write" << endl;
        H5Sclose(memSpaceId);
    }

    dims[0] = totalNumConnectedVessels;
    dims[1] = 1;
    // memSpaceId = H5Screate_simple(2, dims, NULL);

    count[0] = maxNumVesselsPerWrite;

    for (int i = 0; i < numLoopsNeeded; i++)
    {
        // count[0] = subArterialTreeNumberOfVessels / numWritesNeeded;
        memoryOffset = i * ((count[0]) * count[1]) +
            remainingNumberOfVesselsToExport * count[1];
        start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank +
            i * count[0] + remainingNumberOfVesselsToExport;
        //

        memSpaceId = H5Screate_simple(2, count, NULL);
        // cout << "Rank" << mpiRank << ": " << start[0] << ", " << count[0] <<
        // endl;
        H5Sselect_hyperslab(dataSpaceId, H5S_SELECT_SET, start, blockStride,
            count, blockStride);

        H5Dwrite(connectionsArrayDataset, H5T_NATIVE_LLONG, memSpaceId,
            dataSpaceId, dtpl, connectedVesselArray + memoryOffset);
        H5Sclose(memSpaceId);
    }

    // H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);
    H5Pclose(dtpl);
}

void WriteSubConnectedVesselDataToFile(long numVesselsExport,
    long long rootNumVessels, long long totalNumConnectedVessels,
    hid_t numConnectionsDataset, std::vector<short>& numConnectedVesselsArray,
    hid_t connectionsArrayDataset, std::vector<long long>& connectedVesselArray)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    // hid_t dtpl = CreateCollectiveDataTransferPropertiesList();
    // hid_t memSpaceId, dataSpaceId;

    hsize_t dims[2];
    hsize_t start[2], count[2], blockStride[2];
    dims[0] = numVesselsExport;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = rootNumVessels + mpiRank * numVesselsExport;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 1;

    CollectiveHdfBlockExport(numConnectionsDataset, start[0], count[0], numConnectedVesselsArray);
    numConnectedVesselsArray.clear();

    // H5Dclose(numConnectionsDataset);
    // dtpl = CreateCollectiveDataTransferPropertiesList();
    // const long maxNumVesselsPerWrite = 67108864;
    // cout << subArterialTreeNumberOfVessels << ", " << arteries.size() << ", "
    // << maxNumVesselsPerWrite<< endl;

    const long long remainingNumberOfVesselsToExport =
        totalNumConnectedVessels % (4 * maxNumVesselsPerWrite);
    const int numLoopsNeeded =
        totalNumConnectedVessels / (4 * maxNumVesselsPerWrite);
    size_t memoryOffset;

    dims[0] = totalNumConnectedVessels;
    dims[1] = 1;

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank;
    start[1] = 0;
    count[0] = totalNumConnectedVessels;
    count[1] = 1;
    // dataSpaceId = H5Dget_space(connectionsArrayDataset);
    if (remainingNumberOfVesselsToExport > 0)
    {
        count[0] = remainingNumberOfVesselsToExport;

        start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank;
        // start[1] = 0;
        CollectiveHdfBlockExport(connectionsArrayDataset, start[0], count[0], connectedVesselArray);
    }

    dims[0] = totalNumConnectedVessels;
    dims[1] = 1;
    // memSpaceId = H5Screate_simple(2, dims, NULL);

    count[0] = (4 * maxNumVesselsPerWrite);

    for (int i = 0; i < numLoopsNeeded; i++)
    {
        // count[0] = subArterialTreeNumberOfVessels / numWritesNeeded;
        memoryOffset = i * ((count[0]) * count[1]) +
            remainingNumberOfVesselsToExport * count[1];
        start[0] = rootNumVessels * 4 - 4 + totalNumConnectedVessels * mpiRank +
            i * count[0] + remainingNumberOfVesselsToExport;
        //

        CollectiveHdfBlockExport(connectionsArrayDataset, start[0], count[0],
            connectedVesselArray, memoryOffset);
    }
    // H5Sclose(memSpaceId);
    // H5Sclose(dataSpaceId);
    // H5Pclose(dtpl);
}

void ExportSubVesselsToHdf5(hid_t file_id, int totalLevels, int subLevels,
    VesselVector& arteries,
    VesselVector& veins)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    hid_t datasetId, geomDatasetId, nodeDatasetId, dataspaceId, dcpl, 
        conductanceDatasetId, numConnectionsDataset, connectionsArrayDataset;
    herr_t status;

    // cout << "Length of Arteries and Veins" << arteries.size() << " " <<
    // veins.size() << endl;

    long long totalNumVessels = pow(2, totalLevels + 1) - 2;
    long long subNumVessels = pow(2, subLevels + 1) - 1;
    long long rootNumVessels = pow(2, totalLevels - subLevels - 1) - 1;

    // cout << subNumVessels << endl;
    // cout << rootNumVessels << endl;

    double time;

    // cout << "Write Vessel Data" << endl;
    //!time = -omp_get_wtime();
    // if (mpiRank == 0)
    //     cout << "Write Geom Data to File " << endl;

    WriteSubVesselsToFile(subNumVessels, rootNumVessels, file_id, arteries, veins);

    //! if (mpiRank == 0)
    //!     cout << "Write Geom Data time: " << omp_get_wtime() + time << endl;

    vector<short> numConnectedVesselsArray(2 * subNumVessels);
    long long totalNumConnectedVessels = PopulateNumConnectedVesselsArray(arteries, veins, numConnectedVesselsArray, subNumVessels);
    vector<long long> connectedVesselArray(totalNumConnectedVessels);
    PopulateConnectedVesselArray(subNumVessels, arteries, veins, connectedVesselArray);
    numConnectionsDataset = OpenHdfDataset(file_id, GEOM_GROUP_NAME, NUM_CONNECTED_VESSELS_DATASET);
    connectionsArrayDataset = OpenHdfDataset(file_id, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);
    WriteSubConnectedVesselDataToFile(2 * subNumVessels, 2 * rootNumVessels,
        totalNumConnectedVessels, numConnectionsDataset,
        numConnectedVesselsArray, connectionsArrayDataset,
        connectedVesselArray);
    status = H5Dclose(numConnectionsDataset);
    status = H5Dclose(connectionsArrayDataset);
    // delete[] connectedVesselArray;
    connectedVesselArray.clear();
    // delete[] numConnectedVesselsArray;
    numConnectedVesselsArray.clear();
}

void SubVesselExporter(std::string filename,
    VesselVector& arteries,
    VesselVector& veins, int levels, int subLevels)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    // cout << " Rank 0 exporting data" << endl;

    hid_t file_id, group_id, datasetId, geomDatasetId, nodeDatasetId,conductanceDatasetId, dataspaceId, dcpl;
    herr_t status;

    // cout << mpiRank << ": Open File" << endl;
    //! double time = -omp_get_wtime();
    // if (mpiRank == 0)
    //     cout << "Opening file" << endl;

    file_id = OpenHdfFile(filename);

    //! if (mpiRank == 0)
    //!     cout << "time: " << omp_get_wtime() + time << endl;
    // cout << mpiRank << ": File Open " << endl;

    ExportSubVesselsToHdf5(file_id, levels, subLevels, arteries, veins);

    status = H5Fclose(file_id);

    return;
}

void SubVesselExporter(hid_t file_id,
    VesselVector& arteries,
    VesselVector& veins, int levels, int subLevels)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    // cout << " Rank 0 exporting data" << endl;

    //hid_t group_id, datasetId, geomDatasetId, nodeDatasetId, conductanceDatasetId, dataspaceId, dcpl;
    herr_t status;
    
    // cout << mpiRank << ": Open File" << endl;
    //! double time = -omp_get_wtime();
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (mpiRank == 0) cout << "Opening file" << endl;
    // file_id = OpenHdfFile(filename);
    // if (mpiRank == 0) cout << "time: " << omp_get_wtime() + time << endl;
    // cout << mpiRank << ": File Open " << endl;
    
    ExportSubVesselsToHdf5(file_id, levels, subLevels, arteries, veins);

    //! time = -omp_get_wtime();
    // if (mpiRank == 0) cout << "Closing file" << endl;

    status = H5Fclose(file_id);

    // if (mpiRank == 0) cout << "time: " << omp_get_wtime() + time << endl;

    return;
}