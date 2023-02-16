#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <hdf5.h>

#include "GeneratorSharedHdfFunctions.h"
#include "vessel.h"
#include "VesselNetworkExporter.h"
#include "coordinate.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef vector<Vessel> VesselVector;
typedef unsigned long ulong;

#define NUM_ITEMS_PER_CHUNK 1024

void WriteVesselDataToFileSerial(hid_t geomDatasetId, double* geomArray,
    hid_t nodeDatasetId, long long* nodeArray, hid_t conductanceDatasetId, double* conductanceArray)
{
    H5Dwrite(geomDatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        geomArray);
    H5Dwrite(nodeDatasetId, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        nodeArray);
    H5Dwrite(conductanceDatasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
        conductanceArray);
}

void WriteConnectedVesselDataToFileSerial(hid_t numConnectionsDataset,
    short* numConnectedVesselsArray, hid_t connectionsArrayDataset,
    long long* connectedVesselArray)
{
    H5Dwrite(numConnectionsDataset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, numConnectedVesselsArray);

    H5Dwrite(connectionsArrayDataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, connectedVesselArray);

    bool value = true;
    hsize_t dims = 1;
    hid_t attributeDataSpce = H5Screate_simple(1, &dims, NULL);
    hid_t attributeId = H5Acreate(connectionsArrayDataset, "ValidOrdering",
        H5T_NATIVE_HBOOL, attributeDataSpce, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attributeId, H5T_NATIVE_HBOOL, &value);
    H5Sclose(attributeDataSpce);
    H5Aclose(attributeId);
}

void ConvertSequentialVesselToGeometryRow(Vessel& vessel, double* geomArray, long long& geomTableIndex) 
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

void ConvertSequentialVesselToNodeRow(Vessel& vessel, long long* nodeArray,long long& nodeTableIndex)
{
    nodeArray[nodeTableIndex] = vessel.get_nodes()[0];
    nodeArray[nodeTableIndex + 1] = vessel.get_nodes()[1];
}

void ConvertSequentialVesselToConductanceRow(Vessel& vessel, double* conductanceArray,long long& rowIndex)
{
    conductanceArray[rowIndex] = vessel.get_conductance();
}

// Function to serialize vessel data for root vessels in sequential algorithm
void ConvertSequentialVesselVectorsToDataTables(long long numVesselsToWrite, VesselVector& vessels, 
			double* geomArray, long long* nodeArray, double* conductanceArray) 
{
	long long rowIndex = 0, geomTableIndex = 0, nodeTableIndex = 0;
    VesselVector::size_type i;
	
    for (i = 0; i < numVesselsToWrite; i++)
    {
		ConvertSequentialVesselToGeometryRow(vessels[i], geomArray, geomTableIndex);
        ConvertSequentialVesselToNodeRow(vessels[i], nodeArray, nodeTableIndex);
        ConvertSequentialVesselToConductanceRow(vessels[i], conductanceArray, rowIndex);
		geomTableIndex += 7;
        nodeTableIndex += 2;
        rowIndex += 1;
    }
}

void ExportVesselsToHdf5(hid_t file_id, int levels,
    VesselVector& vessels)
{
    hid_t datasetId, geomDatasetId, nodeDatasetId, dataspaceId, dcpl,
        numConnectionsDataset, connectionsArrayDataset, conductanceDatasetId;
    herr_t status;

    long long numVessels = vessels.size();
   
    geomDatasetId = CreateGeometryHdfDataset(file_id, levels, numVessels);
    nodeDatasetId = CreateNodeHdfDataset(file_id, levels, numVessels);
    conductanceDatasetId = CreateHealthyConductanceHdfDataset(file_id, levels, numVessels);

    double* geomArray = new double[numVessels * 7];
    long long* nodeArray = new long long[numVessels * 2];
    double* conductanceArray = new double[numVessels];

    //ConvertVesselVectorsToDataTables(arteries.size(), arteries, veins, geomArray, nodeArray);
    ConvertSequentialVesselVectorsToDataTables(numVessels, vessels, geomArray, nodeArray, conductanceArray);

    WriteVesselDataToFileSerial(geomDatasetId, geomArray, nodeDatasetId, nodeArray, conductanceDatasetId, conductanceArray);

    delete[] nodeArray;
    delete[] geomArray;
    delete[] conductanceArray;

    status = H5Dclose(geomDatasetId);
    status = H5Dclose(nodeDatasetId);
    status = H5Dclose(conductanceDatasetId);

    //!short* numConnectedVesselsArray = new short[numVessels];

    //!long long totalNumConnectedVessels = PopulateNumConnectedVesselsArray(arteries, veins, numConnectedVesselsArray, arteries.size());

    // cout << totalNumConnectedVessels << endl;
    //!long long* connectedVesselArray = new long long[totalNumConnectedVessels];
    // cout << "entering populate connected vessels" << endl;
    //!PopulateConnectedVesselArray(arteries.size(), arteries, veins, connectedVesselArray);

    //!numConnectionsDataset = createNumConnectedVesselsDataset(file_id, levels, numVessels);
    //!connectionsArrayDataset = createConnectedVesselsDataset(file_id, levels, totalNumConnectedVessels);
    // cout << "entering connected data writer" << endl;
    //!WriteConnectedVesselDataToFileSerial(numConnectionsDataset, numConnectedVesselsArray, connectionsArrayDataset, connectedVesselArray);
    // cout << "leaving connected data writer" << endl;

    //! status = H5Dclose(numConnectionsDataset);
    //! status = H5Dclose(connectionsArrayDataset);

    //! delete[] connectedVesselArray;
    //! delete[] numConnectedVesselsArray;
}

void VesselNetworkExporter_HDF5(std::string filename, int levels, VesselVector& vessels, std::vector<double>& boundingBox)
{
    hid_t file_id, group_id, datasetId, geomDatasetId, nodeDatasetId,
        dataspaceId, dcpl;
    herr_t status;

    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    CreateHdfGroups(file_id);

    CreateAndExportNumberOfLevels(file_id, levels);

    CreateAndExportBoundingBox(file_id, boundingBox);

    ExportVesselsToHdf5(file_id, levels, vessels);

    CreateAndExportBoundaryConditions(file_id, levels);

    CreatePartitionAttributes(file_id);

    status = H5Fclose(file_id);

    return;
}