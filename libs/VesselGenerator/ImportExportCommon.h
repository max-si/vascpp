#pragma once

#include <string>

#include <hdf5.h>

#include "Hdf5FileStructureNames.h"
#include "ImportExportHdfBlock.h"

#ifdef CLUSTER
const long long maxNumVesselsPerWrite = (9LL * 67108864LL) / 16LL;
#else
const long long maxNumVesselsPerWrite = 1024 * 1024;
//const  long long maxNumVesselsPerWrite = 4;
#endif

//* Open or Create an HDF5 file
hid_t OpenHdfFile(std::string& fileName);
hid_t OpenHdfFileSerial(std::string& fileName);
hid_t CreateHdfFile(std::string filename);
hid_t CreateHdfFileSerial(std::string filename);
hid_t OpenHdfFileReadOnly(std::string& fileName);
hid_t OpenHdfFileReadOnlySerial(std::string& fileName);

//Open an HDF5 Dataset
hid_t OpenHdfDataset(
    hid_t file_id, std::string groupName, std::string datasetName);

//Import Functions
int RetrieveNumberOfLevelsFromFile(hid_t file_id);
long long ImportConnectedVesselData(
    hid_t fileId, long long*& connectedVesselData);
int GetNumPartitionsInFile(hid_t fileId);
long long GetNumRowsInDataset(hid_t datasetId);

void ImportGeomDataBlock(hid_t datasetId, long long startVesselIndex,
    long long numVesselsToImport, std::vector<double>& data);
void ImportNodeDataBlock(hid_t datasetId, long long startVesselIndex,
    long long numVesselsToImport, std::vector<long long>& data);
void ImportGeometryAndNodeDataBlock(hid_t geomDatasetId, hid_t nodeDatasetId,
    long long startVesselIndex, long long numVesselsToImport,
    std::vector<double>& geomData, std::vector<long long>& nodeData);
void NullImportGeometry(hid_t datasetId);
void NullImportNode(hid_t datasetId);
void PerformNullImportOfGeometryAndNodes(
    hid_t geomDatasetId, hid_t nodeDatasetId);
std::vector<long long> ImportConnectedVesselData(hid_t fileId);

std::vector<long long> ImportNumberOfVesselsPerPartitionData(
    hid_t fileId, int numPartitions);
long long GetTotalNumberOfVessels(hid_t fileId);
//Export Functions
void ExportConnectedVesselData(hid_t fileId, long long*& connectedVesselData);
void ExportConnectedVesselData(
    hid_t fileId, std::vector<long long>& connectedVesselData);

int RetrieveNumberOfLevelsFromFileSerial(hid_t file_id);
//Other Utilities
hid_t CreateCollectiveDataTransferPropertiesList();
hid_t CreateIndependentDataTransferPropertiesList();
hid_t CreateFAPL();
hid_t CreateSerialFAPL();
void DetermineIOLoopCount(long long& remainingNumberOfVesselsToExport,
    int& numLoopsNeeded, int& globalNumLoopsNeeded);