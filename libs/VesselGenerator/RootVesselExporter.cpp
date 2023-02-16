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

using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef std::vector<Vessel> VesselVector;
typedef unsigned long ulong;

#define NUM_ITEMS_PER_CHUNK 1024

void WriteRootVesselsToFile(int numVesselsExport, hid_t geomDatasetId,
    double* geomArray, hid_t nodeDatasetId, long long* nodeArray, hid_t conductanceDatasetId, double* conductanceArray)
{
    hid_t memSpaceId, dataSpaceId;
    hsize_t dims[2];
    hsize_t start[2], count[2], blockStride[2];
    dims[0] = numVesselsExport;
    dims[1] = 7;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(geomDatasetId);

    //write root vessels to file
    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = 0;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 7;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(geomDatasetId, H5T_NATIVE_DOUBLE, memSpaceId, dataSpaceId,
        H5P_DEFAULT, geomArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);

    // write root nodes to file
    dims[0] = numVesselsExport;
    dims[1] = 2;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(nodeDatasetId);

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = 0;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 2;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(nodeDatasetId, H5T_NATIVE_LLONG, memSpaceId, dataSpaceId,
        H5P_DEFAULT, nodeArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);

    // write root conductances to file
    dims[0] = numVesselsExport;
    dims[1] = 1;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(conductanceDatasetId);

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = 0;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 1;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(conductanceDatasetId, H5T_NATIVE_DOUBLE, memSpaceId, dataSpaceId, H5P_DEFAULT, conductanceArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);
}

void WriteRootConnectedVesselDataToFile(int numVesselsExport,
    long long totalNumConnectedVessels, hid_t numConnectionsDataset,
    short* numConnectedVesselsArray, hid_t connectionsArrayDataset,
    long long* connectedVesselArray)
{
    hid_t memSpaceId, dataSpaceId;
    hsize_t dims[2];
    hsize_t start[2], count[2], blockStride[2];
    dims[0] = numVesselsExport;
    dims[1] = 1;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(numConnectionsDataset);

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = 0;
    start[1] = 0;
    count[0] = numVesselsExport;
    count[1] = 1;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(numConnectionsDataset, H5T_NATIVE_SHORT, memSpaceId, dataSpaceId,
        H5P_DEFAULT, numConnectedVesselsArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);

    dims[0] = totalNumConnectedVessels;
    dims[1] = 1;
    memSpaceId = H5Screate_simple(2, dims, NULL);
    dataSpaceId = H5Dget_space(numConnectionsDataset);

    blockStride[0] = 1;
    blockStride[1] = 1;
    start[0] = 0;
    start[1] = 0;
    count[0] = totalNumConnectedVessels;
    count[1] = 1;
    H5Sselect_hyperslab(
        dataSpaceId, H5S_SELECT_SET, start, blockStride, count, blockStride);

    H5Dwrite(connectionsArrayDataset, H5T_NATIVE_LLONG, memSpaceId, dataSpaceId,
        H5P_DEFAULT, connectedVesselArray);
    H5Sclose(memSpaceId);
    H5Sclose(dataSpaceId);

    bool value = true;
    dims[0] = 1;
    hid_t attributeDataSpce = H5Screate_simple(1, dims, NULL);
    hid_t attributeId = H5Acreate(connectionsArrayDataset, "ValidOrdering",
        H5T_NATIVE_HBOOL, attributeDataSpce, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attributeId, H5T_NATIVE_HBOOL, &value);
    H5Sclose(attributeDataSpce);
    H5Aclose(attributeId);
}

void ExportRootVesselsToHdf5(hid_t file_id, int totalLevels, int rootLevels,
    VesselVector& arteries,
    VesselVector& veins)
{
    hid_t datasetId, geomDatasetId, nodeDatasetId, dataspaceId, dcpl,
        numConnectionsDataset, connectionsArrayDataset, conductanceDatasetId;
    herr_t status;

    long long totalNumVessels = pow(2, totalLevels + 1) - 2;
    long long rootNumVessels = pow(2, rootLevels - 1) - 1;
    //cout << "Create Geom Dataset" << endl;
    geomDatasetId = CreateGeometryHdfDataset(file_id, totalLevels, totalNumVessels);
    
    //cout << "Create Node Dataset" << endl;
    nodeDatasetId = CreateNodeHdfDataset(file_id, totalLevels, totalNumVessels);

    //cout << "Create Conductance Dataset" << endl;
    conductanceDatasetId = CreateHealthyConductanceHdfDataset(file_id, totalLevels, totalNumVessels);

    double* geomArray = new double[2 * rootNumVessels * 7];
    long long* nodeArray = new long long[2 * rootNumVessels * 2];
    double* conductanceArray = new double[2 * rootNumVessels];

    // cout << "Populate arrays" << endl;
    ConvertRootVesselVectorsToDataTables(rootNumVessels, arteries, veins, geomArray, nodeArray, conductanceArray);

    // for(int i = 0; i < 2*rootNumVessels; i++) {
    //     std::cout << conductanceArray[i] << " ";
    // } std::cout << std::endl;

    // cout << "Write Vessel Data" << endl;
    WriteRootVesselsToFile(2 * rootNumVessels, geomDatasetId, geomArray, nodeDatasetId, nodeArray, conductanceDatasetId, conductanceArray);

    delete[] nodeArray;
    delete[] geomArray;
    delete[] conductanceArray;

    status = H5Dclose(geomDatasetId);
    status = H5Dclose(nodeDatasetId);
    status = H5Dclose(conductanceDatasetId);

    short* numConnectedVesselsArray = new short[2 * rootNumVessels];
    long long totalNumConnectedVessels = PopulateRootNumConnectedVesselsArray(arteries, veins, numConnectedVesselsArray, rootNumVessels);
    long long* connectedVesselArray = new long long[totalNumConnectedVessels];
    PopulateRootConnectedVesselArray(rootNumVessels, arteries, veins, connectedVesselArray);
    numConnectionsDataset = createNumConnectedVesselsDataset(file_id, totalLevels, totalNumVessels);
    connectionsArrayDataset = createConnectedVesselsDataset(file_id, totalLevels, totalNumConnectedVessels);
    WriteRootConnectedVesselDataToFile(2 * rootNumVessels,
        totalNumConnectedVessels, numConnectionsDataset,
        numConnectedVesselsArray, connectionsArrayDataset,
        connectedVesselArray);
    status = H5Dclose(numConnectionsDataset);
    status = H5Dclose(connectionsArrayDataset);
    delete[] connectedVesselArray;
    delete[] numConnectedVesselsArray;
}

hid_t RootVesselExporter(std::string filename,
    VesselVector& arteries,
    VesselVector& veins, int levels, int rootLevels,
    std::vector<double>& boundingBox)
{
    int mpiSize = 1;
    int mpiRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    hid_t file_id;
    herr_t status;

    // cout << "Create File" << endl;
    file_id = CreateHdfFile(filename);

    // cout << "Create Groups" << endl;
    CreateHdfGroups(file_id);

    // cout << "Export numLevels" << endl;
    CreateAndExportNumberOfLevels(file_id, levels);

    // cout << "Export Boundary Conditions" << endl;
    CreateAndExportBoundaryConditions(file_id, levels);

    // cout << "Export Attributes" << endl;
    CreatePartitionAttributes(file_id);

    // cout << "Export Bounding Box" << endl;
    CreateAndExportBoundingBox(file_id, boundingBox);

    // cout << "Export Vessel Geom & Nodes" << endl;
    ExportRootVesselsToHdf5(file_id, levels, rootLevels, arteries, veins);

    H5Fflush(file_id, H5F_SCOPE_GLOBAL);

    return file_id;
}
