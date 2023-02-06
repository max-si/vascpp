#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>

#include "HdfTypeTemplate.h"

hid_t CreateCollectiveDataTransferPropertiesList();
hid_t CreateIndependentDataTransferPropertiesList();
namespace {

    inline void CheckIfIoWasCollective(hid_t dtplId)
    {
        uint32_t localFailureReason;
        uint32_t globalFailureReason;
        H5Pget_mpio_no_collective_cause(
            dtplId, &localFailureReason, &globalFailureReason);
        if (globalFailureReason == H5D_MPIO_COLLECTIVE ||
            globalFailureReason == H5D_MPIO_SET_INDEPENDENT)
            return;
        if (globalFailureReason & H5D_MPIO_DATATYPE_CONVERSION)
        {
            std::cout << "Collective I/O was not performed because datatype "
                         "conversions were required."
                      << std::endl;
            ;
        }
        if (globalFailureReason & H5D_MPIO_DATA_TRANSFORMS)
        {
            std::cout << "Collective I/O was not performed because data "
                         "transforms needed to be applied."
                      << std::endl;
            ;
        }
        if (globalFailureReason & 8)
        {
            std::cout << "Collective I/O was not performed because the "
                         "selected file driver was MPI-POSIX."
                      << std::endl;
            ;
        }
        if (globalFailureReason & H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES)
        {
            std::cout << "Collective I/O was not performed because one of the "
                         "dataspaces was neither simple nor scalar."
                      << std::endl;
            ;
        }
        if (globalFailureReason & 32)
        {
            std::cout << "Collective I/O was not performed because there were "
                         "point selections in one of the dataspaces."
                      << std::endl;
            ;
        }
        if (globalFailureReason & H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET)
        {
            std::cout << "	Collective I/O was not performed because the "
                         "dataset was neither contiguous nor chunked."
                      << std::endl;
            ;
        }
        //		if (globalFailureReason & H5D_MPIO_FILTERS)
        //		{
        //			std::cout << "Collective I/O was not performed because filters needed to be applied." << std::endl;;
        //		}

        return;
    }
}    // namespace
/**
 * Template based block import of data from HDF5 dataset
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param dataOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void CollectiveHdfBlockImport(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, hsize_t dataOffset = 0)
{
    data.resize(numRowsToImport * numColumns);
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = numRowsToImport;
    count[1] = numColumns;
    start[0] = startRowIndex;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
        count, strideAndBlocks);
    H5Dread(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        data.data() + dataOffset);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Template based null import, used for collective imports with no data
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 */
template <typename T, int dims = 1, int numColumns = 1>
void CollectiveHdfNullImport(hid_t datasetId)
{
    double* nullDataPointer = nullptr;
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = 0;
    count[1] = numColumns;
    start[0] = 0;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    hid_t dataspaceId = H5Dget_space(datasetId);
    H5Sselect_none(memspaceId);
    H5Sselect_none(dataspaceId);

    H5Dread(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        nullDataPointer);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Indepndent import using a template based format.
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param dataOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void IndependentHdfBlockImport(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, hsize_t dataOffset = 0)
{
    data.resize(numRowsToImport * numColumns);
    hid_t dtplId = CreateIndependentDataTransferPropertiesList();

    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = numRowsToImport;
    count[1] = numColumns;
    start[0] = startRowIndex;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
        count, strideAndBlocks);
    H5Dread(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        data.data() + dataOffset);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * A speciallized template importer for small dataset, where MPI_BCAST can be used to import lots of small data
 * @tparam T
 * @tparam mpiType
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param dataOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void HdfBlockImportBcast(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, MPI_Datatype mpiType, hsize_t dataOffset = 0)
{
    data.resize(numRowsToImport * numColumns);
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0)
    {
        hid_t dtplId = CreateIndependentDataTransferPropertiesList();
        hid_t dataspaceId = H5Dget_space(datasetId);
        hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
        count[0] = numRowsToImport;
        count[1] = numColumns;
        start[0] = startRowIndex;
        start[1] = 0;
        hid_t memspaceId = H5Screate_simple(dims, count, NULL);
        H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
            count, strideAndBlocks);
        H5Dread(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
            data.data() + dataOffset);
        CheckIfIoWasCollective(dtplId);
        H5Sclose(dataspaceId);
        H5Sclose(memspaceId);
        H5Pclose(dtplId);
    }

    MPI_Bcast(
        data.data(), numColumns * numRowsToImport, mpiType, 0, MPI_COMM_WORLD);
}
// Export functions

/**
 * Template based collective block export of HDF5 Data
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param dataOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void CollectiveHdfBlockExport(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, hsize_t dataOffset = 0)
{
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();
    //std::cout << numRowsToImport << " " << startRowIndex <<" "<< data.size() << std::endl;
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = numRowsToImport;
    count[1] = numColumns;
    start[0] = startRowIndex;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
        count, strideAndBlocks);

    H5Dwrite(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        data.data() + dataOffset);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Template based Export for when no Data needs to be exported by a compute node
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 */
template <typename T, int dims = 1, int numColumns = 1>
void CollectiveHdfNullExport(hid_t datasetId)
{
    double* nullDataPointer = nullptr;
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = 0;
    count[1] = numColumns;
    start[0] = 0;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    hid_t dataspaceId = H5Dget_space(datasetId);
    H5Sselect_none(memspaceId);
    H5Sselect_none(dataspaceId);
    H5Dwrite(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        nullDataPointer);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Template-based input for independent inport of a parallel file
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param dataOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void IndependentHdfBlockExport(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, hsize_t dataOffset = 0)
{
    hid_t dtplId = CreateIndependentDataTransferPropertiesList();

    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = numRowsToImport;
    count[1] = numColumns;
    start[0] = startRowIndex;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
        count, strideAndBlocks);
    H5Dwrite(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        data.data() + dataOffset);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Import template function that imports data at the end of a data vector
 * @tparam T
 * @tparam dims
 * @tparam numColumns
 * @param datasetId
 * @param startRowIndex
 * @param numRowsToImport
 * @param data
 * @param vectorOffset
 */
template <typename T, int dims = 1, int numColumns = 1>
void CollectiveHdfAppendingImport(hid_t datasetId, hsize_t startRowIndex,
    hsize_t numRowsToImport, std::vector<T>& data, hsize_t vectorOffset = 0)
{
    //data.resize(numRowsToImport * numColumns + 10);
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = numRowsToImport;
    count[1] = numColumns;
    start[0] = startRowIndex;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(dims, count, NULL);
    H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
        count, strideAndBlocks);
    H5Dread(datasetId, HdfType<T>(), memspaceId, dataspaceId, dtplId,
        data.data() + vectorOffset);
    CheckIfIoWasCollective(dtplId);
    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

/**
 * Specialized import function that used MPI_Broadcast for small datasets
 * @tparam T
 * @tparam mpiType
 * @tparam dims
 * @tparam numColumns
 * @param fileId
 * @param groupName
 * @param datasetName
 * @return
 */
template <typename T, int dims = 1, int numColumns = 1>
std::vector<T> ImportAndBcastSmallDataset(
    hid_t fileId, std::string groupName, std::string datasetName, MPI_Datatype mpiType)
{
    int mpiRank, mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    std::string datasetPath = groupName + datasetName;

    hid_t datasetId = H5Dopen2(fileId, datasetPath.c_str(), H5P_DEFAULT);
    hid_t dataspaceId = H5Dget_space(datasetId);
    hsize_t numBcs = 0;
    H5Sget_simple_extent_dims(dataspaceId, &numBcs, NULL);
    std::vector<T> bcNodes(numBcs);
    H5Sclose(dataspaceId);

    if (!mpiRank)
    {
        IndependentHdfBlockImport<T, dims>(datasetId, 0, numBcs, bcNodes);
    }

    MPI_Bcast(bcNodes.data(), numBcs, mpiType, 0, MPI_COMM_WORLD);

    H5Dclose(datasetId);

    return std::move(bcNodes);
};