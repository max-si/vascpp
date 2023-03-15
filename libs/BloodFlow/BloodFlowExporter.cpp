#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "BloodFlowVessel.h"
#include "NetworkDescription.h"
//#include <vascular/core/GlobalDefs.h>
#include "VesselGenerator/ImportExportCommon.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef vector<BloodFlowVessel> VesselVector;
typedef unsigned long ulong;

void ExportFlowsToFile(const NetworkDescription& network, hid_t datasetId, double (BloodFlowVessel::*flowTypePointer)() const);

hid_t CreateHdfFlowDataset(hid_t fileId, const NetworkDescription& network, std::string groupName, std::string datasetName);

void BloodFlowExporter(hid_t fileId, const NetworkDescription& network)
{
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    //export healthy flows
    hid_t datasetId = CreateHdfFlowDataset(fileId, network, FLOW_GROUP_NAME, HEALTHY_FLOW_DATASET);
    double (BloodFlowVessel::*flowTypePointer)() const = &BloodFlowVessel::GetHealthyFlow;
    ExportFlowsToFile(network, datasetId, flowTypePointer);
    H5Dclose(datasetId);

    //!export damaged flows
    // datasetId = CreateHdfFlowDataset(fileId, network, FLOW_GROUP_NAME, DAMAGED_FLOW_DATASET);
    // flowTypePointer = &BloodFlowVessel::GetDamagedFlow;
    // ExportFlowsToFile(network, datasetId, flowTypePointer);
    // H5Dclose(datasetId);

    H5Fflush(fileId, H5F_SCOPE_GLOBAL);
}

/**
 *  \breif Creates blood flow dataset
 *
 *  Function creates a dataset for blood flow rates in the HDF5 file
 *  Could eventually be replaced with the CreateHDFDataset functions
 *
 * @param[in] fileId HDF5 file handle to be exported to
 * @param[in] network Network Definition (only used for size of network)
 * @param[in] groupName Name of group to contain dataset
 * @param[in] datasetName Name of Dataset
 * @return
 */
hid_t CreateHdfFlowDataset(hid_t fileId, const NetworkDescription& network, std::string groupName, std::string datasetName)
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
    if (datasetId > 0) {
        return datasetId;
    }

    dims[0] = network.getTotalNumberOfVessels();
    dims[1] = 1;
    dataspaceId = H5Screate_simple(1, dims, NULL);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    datasetId = H5Dcreate2(fileId, datasetPath.c_str(), H5T_NATIVE_DOUBLE, dataspaceId, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Pclose(dcpl);
    status = H5Sclose(dataspaceId);

    return datasetId;
}

/**
 *  \breif Function that exports data to file
 * @param[in] network Vascular Network to be exported
 * @param[in] datasetId HDF5 file handle for Dataset
 * @param[in] flowTypePointer Pointer to method to extract data from class
 */
void ExportFlowsToFile(const NetworkDescription& network, hid_t datasetId,
    double (BloodFlowVessel::*flowTypePointer)() const)    //TODO Update to Block I/O
{
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    //Serialize Data from vessel objects
    vector<double> dataArray(network.getLocalNumberOfVessels());
    const vector<BloodFlowVessel>& vessels = network.getLocalVesselVector();
    ulong arrayIndex = 0;

    VesselVector::const_iterator itr;
    for (itr = vessels.begin(); itr != vessels.end(); itr++)
        dataArray[arrayIndex++] = ((*itr).*flowTypePointer)();

    //Write data to file
    hid_t dataTransferPropertiesList = CreateCollectiveDataTransferPropertiesList();

    hsize_t count = vessels.size();
    hsize_t offset = network.getGlobalVesselStartIndex();

    //cout << "Rank " << mpiRank << " offset: " << offset << " count: " << count << " numVessels: "<<vessels.size() << endl;

    hid_t memSpace = H5Screate_simple(1, &count, NULL);
    hid_t fileSpace = H5Dget_space(datasetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &offset, NULL, &count, NULL);

    H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, memSpace, fileSpace, dataTransferPropertiesList, dataArray.data());
    H5Pclose(dataTransferPropertiesList);
    //free(dataArray);
    return;
}