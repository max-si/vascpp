#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <hdf5.h>


#include "Hdf5FileStructureNames.h"
#include "GeneratorSharedHdfFunctions.h"
#include "vessel.h"
#include "VesselNetworkExporter.h"
#include "coordinate.h"

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::vector;

typedef std::vector<Vessel> VesselVector;
typedef unsigned long ulong;

#define NUM_ITEMS_PER_CHUNK 512

void ConvertBoundingBoxToArray(
    std::vector<double> boundingBox, double* boundingBoxArray)
{
    boundingBoxArray[0] = boundingBox[0];
    boundingBoxArray[1] = boundingBox[1];
    boundingBoxArray[2] = boundingBox[2];
    boundingBoxArray[3] = boundingBox[3];
    boundingBoxArray[4] = boundingBox[4];
    boundingBoxArray[5] = boundingBox[5];
}

void ConvertVesselToGeometryRow(
    Vessel& vessel, double* geomArray, long long& geomTableIndex)
{
    const Coordinate& start = vessel.get_startingPoint();
    const Coordinate& end = vessel.get_endingPoint();

    geomArray[geomTableIndex] = start.x;
    geomArray[geomTableIndex + 1] = start.y;
    geomArray[geomTableIndex + 2] = start.z;
    geomArray[geomTableIndex + 3] = end.x;
    geomArray[geomTableIndex + 4] = end.y;
    geomArray[geomTableIndex + 5] = end.z;
    geomArray[geomTableIndex + 6] = vessel.get_radius();
}

void ConvertVesselToNodeRow(Vessel& vessel, long long* nodeArray,
    long long& nodeTableIndex)
{
    // nodeArray[nodeTableIndex] = vessel.get_nodes()[0];
    // nodeArray[nodeTableIndex + 1] = vessel.get_nodes()[1];
    nodeArray[nodeTableIndex] = vessel.getNode1();
    nodeArray[nodeTableIndex + 1] = vessel.getNode2();
}

void ConvertVesselToConductanceRow(Vessel& vessel, double* conductanceArray, long long& rowIndex)
{
    conductanceArray[rowIndex] = vessel.get_conductance();
}

//*
void ConvertVesselToTableRow(Vessel& vessel, double* geomArray,
    long long& geomTableIndex, long long* nodeArray, long long& nodeTableIndex)
{
    const Coordinate& start = vessel.get_startingPoint();
    const Coordinate& end = vessel.get_endingPoint();

    geomArray[geomTableIndex++] = start.x;
    geomArray[geomTableIndex++] = start.y;
    geomArray[geomTableIndex++] = start.z;
    geomArray[geomTableIndex++] = end.x;
    geomArray[geomTableIndex++] = end.y;
    geomArray[geomTableIndex++] = end.z;
    geomArray[geomTableIndex++] = vessel.get_radius();

    // nodeArray[nodeTableIndex++] = vessel.get_nodes()[0];
    // nodeArray[nodeTableIndex++] = vessel.get_nodes()[1];
    nodeArray[nodeTableIndex++] = vessel.getNode1();
    nodeArray[nodeTableIndex++] = vessel.getNode2();
} /**/

void ConvertRootVesselVectorsToDataTables(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, double* geomArray,
    long long* nodeArray, double* conductanceArray)
{
    long long rowIndex = 0, geomTableIndex = 0, nodeTableIndex = 0;
    VesselVector::size_type i;

    for (i = 0; i < numVesselsToWrite; i++)
    {
        ConvertVesselToGeometryRow(arteries[i], geomArray, geomTableIndex);
        ConvertVesselToNodeRow(arteries[i], nodeArray, nodeTableIndex);
        ConvertVesselToConductanceRow(arteries[i], conductanceArray, rowIndex);
        geomTableIndex += 7;
        nodeTableIndex += 2;
        rowIndex += 1;
    }
    for (i = numVesselsToWrite + 1; i < veins.size(); i++)
    //!for (i = 0; i < numVesselsToWrite; i++)
    {
        ConvertVesselToGeometryRow(veins[i], geomArray, geomTableIndex);
        ConvertVesselToNodeRow(veins[i], nodeArray, nodeTableIndex);
        ConvertVesselToConductanceRow(veins[i], conductanceArray, rowIndex);
        geomTableIndex += 7;
        nodeTableIndex += 2;
        rowIndex += 1;
    }
}

void ConvertVesselVectorsToDataTables(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, double* geomArray,
    long long* nodeArray)
{
    long long rowIndex = 0, geomTableIndex = 0, nodeTableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++)
    {
        ConvertVesselToTableRow(
            arteries[i], geomArray, geomTableIndex, nodeArray, nodeTableIndex);
    }
    for (i = 0; i < veins.size(); i++)
    {
        ConvertVesselToTableRow(
            veins[i], geomArray, geomTableIndex, nodeArray, nodeTableIndex);
    }
}

void ConvertVesselVectorsToGeometryArray(long long numVesselsToWrite,
    VesselVector& vessels, double* geomArray)
{
    long long rowIndex = 0, geomTableIndex = 0;
    VesselVector::size_type i;
#pragma omp parallel for
    for (i = 0; i < numVesselsToWrite; i++)
    {
        geomTableIndex = i * 7;
        ConvertVesselToGeometryRow(vessels[i], geomArray, geomTableIndex);
    }
}

void ConvertVesselVectorsToNodeArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* nodeArray)
{
    long long rowIndex = 0, nodeTableIndex = 0;
    VesselVector::size_type i;
#pragma omp parallel for
    for (i = 0; i < numVesselsToWrite; i++)
    {
        nodeTableIndex = i * 2;
        ConvertVesselToNodeRow(arteries[i], nodeArray, nodeTableIndex);
    }
#pragma omp parallel for
    for (i = 0; i < numVesselsToWrite; i++)
    {
        nodeTableIndex = (numVesselsToWrite + i) * 2;
        ConvertVesselToNodeRow(veins[i], nodeArray, nodeTableIndex);
    }
}

void ConvertVesselVectorsToConductanceArray(long long numVesselsToWrite,
    VesselVector& arteries, VesselVector& veins, double* conductanceArray)
{
    long long rowIndex = 0;
    VesselVector::size_type i;
    for(i = 0; i < numVesselsToWrite; i++) {
        rowIndex = i;
        ConvertVesselToConductanceRow(arteries[i], conductanceArray, rowIndex);
    }
    for (i = 0; i < numVesselsToWrite; i++) {
        rowIndex = numVesselsToWrite + i;
        ConvertVesselToConductanceRow(veins[i], conductanceArray, rowIndex);
    }
}

void CreateHdfGroups(hid_t file_id)
{
    hid_t group_id;
    group_id = H5Gcreate(file_id, DOSE_GROUP_NAME.c_str(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, DAMAGE_GROUP_NAME.c_str(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, FLOW_GROUP_NAME.c_str(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, GEOM_GROUP_NAME.c_str(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
}

void CreateAndExportNumberOfLevels(hid_t file_id, int levels)
{
    int mpiRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    hid_t dataspaceId, dcpl, datasetId, memSpaceId;
    herr_t status;
    hsize_t dims[2];
    dims[0] = 1;
    dims[1] = 1;

    dataspaceId = H5Screate_simple(1, dims, NULL);
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    // status = H5Pset_layout(dcpl, H5D_COMPACT);

    datasetId = H5Dcreate2(file_id, NUM_LEVELS_DATASET.c_str(), H5T_STD_U32LE,
        dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    if (!mpiRank)
        memSpaceId = H5Screate_simple(1, dims, NULL);
    else
    {
        H5Sclose(dataspaceId);
        dataspaceId = H5Screate(H5S_NULL);
        memSpaceId = H5Screate(H5S_NULL);
    }
    hid_t dtplId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dtplId, H5FD_MPIO_INDEPENDENT);

    status = H5Dwrite(
        datasetId, H5T_NATIVE_INT, memSpaceId, dataspaceId, dtplId, &levels);
    H5Pclose(dtplId);
    H5Sclose(memSpaceId);
    status = H5Pclose(dcpl);
    status = H5Dclose(datasetId);
    status = H5Sclose(dataspaceId);
}

void ExportBoundingBox(hid_t file_id, double* boundingBoxArray)
{
    int mpiRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t dataspaceId, dcpl, datasetId, memSpaceId;
    herr_t status;
    hsize_t dims[2];
    dims[0] = 6;
    dims[1] = 1;
    dataspaceId = H5Screate_simple(1, dims, NULL);
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    // status = H5Pset_layout(dcpl, H5D_COMPACT);
    datasetId = H5Dcreate2(file_id, BOUNDING_BOX_DATASET.c_str(),
        H5T_NATIVE_DOUBLE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (!mpiRank)
        memSpaceId = H5Screate_simple(1, dims, NULL);
    else
    {
        H5Sclose(dataspaceId);
        dataspaceId = H5Screate(H5S_NULL);
        memSpaceId = H5Screate(H5S_NULL);
    }
    hid_t dtplId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dtplId, H5FD_MPIO_INDEPENDENT);

    status = H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, memSpaceId, dataspaceId,
        H5P_DEFAULT, boundingBoxArray);
    H5Pclose(dtplId);
    H5Sclose(memSpaceId);
    status = H5Pclose(dcpl);
    status = H5Dclose(datasetId);
    status = H5Sclose(dataspaceId);
}

void CreateAndExportBoundingBox(hid_t file_id, std::vector<double>& boundingBox)
{
    double boundingBoxArray[6];
    ConvertBoundingBoxToArray(boundingBox, boundingBoxArray);

    ExportBoundingBox(file_id, boundingBoxArray);
}

hid_t CreateGeometryHdfDataset(hid_t file_id, int levels, long long numVessels)
{
    hid_t geomDataspaceId, geomDatasetId;
    herr_t status;
    hid_t dcpl;
    hsize_t dims[2];

    dims[0] = numVessels;
    dims[1] = 7;
    geomDataspaceId = H5Screate_simple(2, dims, NULL);
    string datasetPath = GEOM_GROUP_NAME;
    datasetPath += GEOM_DATASET;

    long long numBytesTotal = (7 * 8) * numVessels;
    long long numBytesPerChunk = numBytesTotal;
    if (numBytesTotal > (1024LL * 1024LL * 1024LL))
    {
        numBytesPerChunk = 1024 * 1024 * 1024 / 16;
    }

    int numChunks = 64;

    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[2];
    if (levels > 18)
    {
        chunk[0] = numVessels / numChunks;
        chunk[1] = 7;
        // status = H5Pset_chunk(dcpl, 2, chunk);
    }
    else if (levels > 4)
    {
        chunk[0] = 16;
        chunk[1] = 7;
        // status = H5Pset_chunk(dcpl, 2, chunk);
    }
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    geomDatasetId = H5Dcreate2(file_id, datasetPath.c_str(), H5T_NATIVE_DOUBLE,
        geomDataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(geomDataspaceId);

    return geomDatasetId;
}

hid_t CreateNodeHdfDataset(hid_t file_id, int levels, long long numVessels)
{
    herr_t status;
    hid_t dcpl;
    hid_t nodeDataspaceId, nodeDatasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;
    dims[0] = numVessels;
    dims[1] = 2;
    nodeDataspaceId = H5Screate_simple(2, dims, NULL);
    datasetPath = GEOM_GROUP_NAME;
    datasetPath += NODE_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    long long numBytesTotal = (2 * 8) * numVessels;
    int numChunks = 64;
    // int numChunks = numBytesTotal / (4LL * 1024L * 1024LL);

    if (levels > 18)
    {
        chunk[0] = numVessels / numChunks;
        chunk[1] = 2;
        // status = H5Pset_chunk(dcpl, 2, chunk);
    }
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    nodeDatasetId = H5Dcreate2(file_id, datasetPath.c_str(), H5T_NATIVE_LLONG,
        nodeDataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(nodeDataspaceId);
    return nodeDatasetId;
}

hid_t CreateHealthyConductanceHdfDataset(hid_t file_id, int levels, long long numVessels)
{
    herr_t status;
    hid_t dcpl;
    hid_t conductanceDataspaceId, conductanceDatasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;
    dims[0] = numVessels;
    dims[1] = 1;
    conductanceDataspaceId = H5Screate_simple(2, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += HEALTHY_CONDUCTANCE_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    long long numBytesTotal = 8 * numVessels;
    int numChunks = 64;
    // int numChunks = numBytesTotal / (4LL * 1024L * 1024LL);

    if (levels > 18)
    {
        chunk[0] = numVessels / numChunks;
        chunk[1] = 2;
        // status = H5Pset_chunk(dcpl, 2, chunk);
    }
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);

    conductanceDatasetId = H5Dcreate2(file_id, datasetPath.c_str(), H5T_NATIVE_DOUBLE,
        conductanceDataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(conductanceDataspaceId);
    return conductanceDatasetId;
}

hid_t createNumConnectedVesselsDataset(
    hid_t file_id, int levels, long long numVessels)
{
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;
    dims[0] = numVessels;
    dims[1] = 1;
    dataspaceId = H5Screate_simple(2, dims, NULL);
    datasetPath = GEOM_GROUP_NAME;
    datasetPath += NUM_CONNECTED_VESSELS_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    long long numBytesTotal = (8) * numVessels;
    // int numChunks = numBytesTotal / (4LL * 1024L * 1024LL);
    int numChunks = 64;
    if (levels > 18)
    {
        chunk[0] = numVessels / numChunks;
        chunk[1] = 1;
        // status = H5Pset_chunk(dcpl, 2, chunk);
    }
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    datasetId = H5Dcreate2(file_id, datasetPath.c_str(), H5T_NATIVE_SHORT,
        dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(dataspaceId);
    return datasetId;
}

hid_t createConnectedVesselsDataset(hid_t file_id, int levels,
    vector<long long>::size_type totalNumConnectedVessels)
{
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId;
    hsize_t dims[2], chunk[2];
    string datasetPath;
    dims[0] = pow(2, levels) * 7 - 12;
    dims[1] = 1;
    dataspaceId = H5Screate_simple(2, dims, NULL);
    datasetPath = GEOM_GROUP_NAME;
    datasetPath += CONNECTED_VESSELS_DATASET;
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    long long numBytesTotal = 8 * totalNumConnectedVessels;
    int numChunks = 64;

    if (levels > 18)
    {
        chunk[0] = dims[0] / numChunks;
        chunk[1] = 1;
        // status = H5Pset_chunk(dcpl, 2, chunk);
        if (chunk[0] > dims[0])
        {
            cout << "Error: Chunk > dim" << endl;
            cout << chunk[0] << " > " << dims[0] << endl;
        }
    }

    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY);
    datasetId = H5Dcreate2(file_id, datasetPath.c_str(), H5T_NATIVE_LLONG,
        dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(dataspaceId);
    return datasetId;
}

long long PopulateRootNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins, short* numConnectedVesselsArray,
    long long numVessels)
{
    long long rowIndex = 0, tableIndex = 0;
    long long totalNumberOfConnectedVessels = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = arteries[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    //!for (i = 0; i < numVessels; i++, tableIndex++)
    for (i = numVessels + 1; i < veins.size(); i++, tableIndex++)
    {
        short numConnectedToThisVessel = veins[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }

    return totalNumberOfConnectedVessels;
}

long long PopulateNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins, short* numConnectedVesselsArray,
    long long numVessels)
{
    long long rowIndex = 0, tableIndex = 0;
    long long totalNumberOfConnectedVessels = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = arteries[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    for (i = 0; i < veins.size(); i++, tableIndex++)
    {
        short numConnectedToThisVessel = veins[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    // cout << tableIndex << endl;
    /*
  VesselVector::const_reverse_iterator veinItr;
  for (veinItr = veins.rbegin(); veinItr != veins.rend();
  veinItr++,tableIndex++) { short numConnectedToThisVessel =
  veinItr->getNumConnectedVessels(); numConnectedVesselsArray[tableIndex] =
  numConnectedToThisVessel; totalNumberOfConnectedVessels +=
  numConnectedToThisVessel;
  }/**/
    // cout << "Num connections array:" << endl;
    // for (int i = 0; i < numVessels; i++)
    // cout << numConnectedVesselsArray[i] << endl;

    return totalNumberOfConnectedVessels;
}

long long PopulateNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins,
    std::vector<short>& numConnectedVesselsArray, long long numVessels)
{
    long long rowIndex = 0, tableIndex = 0;
    long long totalNumberOfConnectedVessels = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVessels; i++, tableIndex++)
    {
        short numConnectedToThisVessel = arteries[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    for (i = 0; i < veins.size(); i++, tableIndex++)
    {
        short numConnectedToThisVessel = veins[i].get_numConnectedVessels();
        numConnectedVesselsArray[tableIndex] = numConnectedToThisVessel;
        totalNumberOfConnectedVessels += numConnectedToThisVessel;
    }
    // cout << tableIndex << endl;
    /*
  VesselVector::const_reverse_iterator veinItr;
  for (veinItr = veins.rbegin(); veinItr != veins.rend();
  veinItr++,tableIndex++) { short numConnectedToThisVessel =
  veinItr->getNumConnectedVessels(); numConnectedVesselsArray[tableIndex] =
  numConnectedToThisVessel; totalNumberOfConnectedVessels +=
  numConnectedToThisVessel;
  }/**/
    // cout << "Num connections array:" << endl;
    // for (int i = 0; i < numVessels; i++)
    // cout << numConnectedVesselsArray[i] << endl;

    return totalNumberOfConnectedVessels;
}

void PopulateRootConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* connectedVesselArray)
{
    unsigned long long tableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++)
    {
        const vector<long long>& connectedVessels =
            arteries[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }

    //!for (i = 0; i < numVesselsToWrite; i++)
    for (i = numVesselsToWrite + 1; i < veins.size(); i++)
    {
        const vector<long long>& connectedVessels =
            veins[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
    // cout << tableIndex << endl;
    /*
  VesselVector::const_reverse_iterator veinItr;
  for (veinItr = veins.rbegin(); veinItr != veins.rend(); veinItr++) {
  const vector<unsigned long>& connectedVessels =
  veinItr->getConnectedVessels(); for (vector<unsigned long>::const_iterator itr
  = connectedVessels.begin(); itr != connectedVessels.end(); itr++)
  {
  connectedVesselArray[tableIndex++] = *itr;
  }
  }/**/
    // cout << tableIndex << endl;
}

void PopulateConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins,
    std::vector<long long>& connectedVesselArray)
{
    unsigned long long tableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++)
    {
        const vector<long long>& connectedVessels =
            arteries[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }

    for (i = 0; i < veins.size(); i++)
    {
        const vector<long long>& connectedVessels =
            veins[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
    // cout << tableIndex << endl;
    /*
  VesselVector::const_reverse_iterator veinItr;
  for (veinItr = veins.rbegin(); veinItr != veins.rend(); veinItr++) {
  const vector<unsigned long>& connectedVessels =
  veinItr->getConnectedVessels(); for (vector<unsigned long>::const_iterator itr
  = connectedVessels.begin(); itr != connectedVessels.end(); itr++)
  {
  connectedVesselArray[tableIndex++] = *itr;
  }
  }/**/
    // cout << tableIndex << endl;
}
void PopulateConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* connectedVesselArray)
{
    unsigned long long tableIndex = 0;
    VesselVector::size_type i;
    for (i = 0; i < numVesselsToWrite; i++)
    {
        const vector<long long>& connectedVessels =
            arteries[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }

    for (i = 0; i < veins.size(); i++)
    {
        const vector<long long>& connectedVessels =
            veins[i].get_connectedVessels();
        for (vector<long long>::const_iterator itr = connectedVessels.begin();
             itr != connectedVessels.end(); itr++)
        {
            connectedVesselArray[tableIndex++] = *itr;
        }
    }
    // cout << tableIndex << endl;
    /*
  VesselVector::const_reverse_iterator veinItr;
  for (veinItr = veins.rbegin(); veinItr != veins.rend(); veinItr++) {
  const vector<unsigned long>& connectedVessels =
  veinItr->getConnectedVessels(); for (vector<unsigned long>::const_iterator itr
  = connectedVessels.begin(); itr != connectedVessels.end(); itr++)
  {
  connectedVesselArray[tableIndex++] = *itr;
  }
  }/**/
    // cout << tableIndex << endl;
}

void CreatePartitionAttributes(hid_t fileId)
{
    // open group
    string attributeName = "Num_Partitions";
    hsize_t dims = 1;

    int numParts = 1;

    hid_t dataSpaceId = H5Screate_simple(1, &dims, NULL);

    hid_t attributeId = H5Acreate2(fileId, attributeName.c_str(),
        H5T_NATIVE_INT, dataSpaceId, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attributeId, H5T_NATIVE_INT, &numParts);

    H5Sclose(dataSpaceId);
    H5Aclose(attributeId);

    attributeName = "Is_Partitioned";
    dataSpaceId = H5Screate_simple(1, &dims, NULL);

    int wasPartitioned = 0;

    attributeId = H5Acreate2(fileId, attributeName.c_str(), H5T_NATIVE_INT,
        dataSpaceId, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attributeId, H5T_NATIVE_INT, &wasPartitioned);

    H5Sclose(dataSpaceId);
    H5Aclose(attributeId);
}

double CalculateBloodFlowBoundaryCondition(int levels)
{
    return pow(2.0f, static_cast<double>(--levels)) * M_PI * 0.0025 * 0.0025;
}

void CreateAndExportBoundaryConditions(hid_t file_id, int levels)
{
    int mpiRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    long long boundaryConditionNodes[2];
    short bcType[2];
    bcType[0] = -1;
    bcType[1] = 1;
    double boundaryConditions[2] = {1, 1};
    boundaryConditions[1] = CalculateBloodFlowBoundaryCondition(levels);
    boundaryConditionNodes[0] = 0;
    boundaryConditionNodes[1] = 3 * std::pow(2, levels - 1) - 1;
    // cout << numNodes << endl;

    // create datasets
    herr_t status;
    hid_t dcpl;
    hid_t dataspaceId, datasetId, memSpaceId;
    hsize_t dims[2], chunk[2];
    string datasetPath;
    dims[0] = 2;
    dims[1] = 1;
    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += BOUNDARY_CONDITION_NODE_DATASET;
    hid_t nodesDatasetId = H5Dcreate2(file_id, datasetPath.c_str(),
        H5T_NATIVE_LLONG, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Sclose(dataspaceId);

    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += BOUNDARY_CONDITION_VALUE_DATASET;
    hid_t valuesDatasetId = H5Dcreate2(file_id, datasetPath.c_str(),
        H5T_NATIVE_DOUBLE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Sclose(dataspaceId);

    dataspaceId = H5Screate_simple(1, dims, NULL);
    datasetPath = FLOW_GROUP_NAME;
    datasetPath += BOUNDARY_CONDITION_TYPE_DATASET;
    hid_t typeDatasetId = H5Dcreate2(file_id, datasetPath.c_str(),
        H5T_NATIVE_SHORT, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (!mpiRank)
        memSpaceId = H5Screate_simple(1, dims, NULL);
    else
    {
        H5Sclose(dataspaceId);
        dataspaceId = H5Screate(H5S_NULL);
        memSpaceId = H5Screate(H5S_NULL);
    }
    hid_t dtplId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dtplId, H5FD_MPIO_INDEPENDENT);

    H5Dwrite(nodesDatasetId, H5T_NATIVE_LLONG, memSpaceId, dataspaceId,
        H5P_DEFAULT, boundaryConditionNodes);
    H5Dwrite(valuesDatasetId, H5T_NATIVE_DOUBLE, memSpaceId, dataspaceId,
        H5P_DEFAULT, boundaryConditions);
    H5Dwrite(typeDatasetId, H5T_NATIVE_SHORT, memSpaceId, dataspaceId,
        H5P_DEFAULT, bcType);

    status = H5Sclose(dataspaceId);
    H5Pclose(dtplId);
    H5Sclose(memSpaceId);
    H5Dclose(nodesDatasetId);
    H5Dclose(valuesDatasetId);
    H5Dclose(typeDatasetId);
}

