#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "ImportExportCommon.h"
#include "ImportExportHdfBlock.h"

using std::cout;
using std::endl;
using std::string;

hid_t CreateCollectiveDataTransferPropertiesList()
{
    hid_t dtpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_COLLECTIVE);
    // H5Pset_dxpl_mpio_chunk_opt(dtpl, H5FD_MPIO_CHUNK_MULTI_IO);
    // H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_INDEPENDENT);
    // H5Pset_buffer(dtpl,	1024*1024*256,NULL, NULL);
    // H5Pset_hyper_vector_size(dtpl,2048);

    return dtpl;
}

hid_t CreateIndependentDataTransferPropertiesList()
{
    hid_t dtpl = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_COLLECTIVE);
    H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_INDEPENDENT);
    // H5Pset_buffer(dtpl, 256 * 1024 * 1024, NULL, NULL);

    return dtpl;
}

hid_t OpenHdfDataset(
    hid_t file_id, std::string groupName, std::string datasetName)
{
    herr_t status;

    string datasetPath = groupName + datasetName;

    hid_t fileId = H5Dopen2(file_id, datasetPath.c_str(), H5P_DEFAULT);

    return fileId;
}

hid_t OpenHdfDatasetParallel(
    hid_t file_id, std::string groupName, std::string datasetName)
{
    herr_t status;

    string datasetPath = groupName + datasetName;

    hid_t dtpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_COLLECTIVE);
    // H5Pset_dxpl_mpio(dtpl, H5FD_MPIO_INDEPENDENT);
    H5Pset_buffer(dtpl, 67108864, NULL, NULL);

    return H5Dopen2(file_id, datasetPath.c_str(), dtpl);
}

hid_t OpenHdfFile(std::string& fileName)
{
    hid_t fileAccessPropertiesList = CreateFAPL();

    hid_t fileId =
        H5Fopen(fileName.c_str(), H5F_ACC_RDWR, fileAccessPropertiesList);

    if (fileId < 0)
    {
        cout << "File " << fileName << " Not found in local directory" << endl;
        throw;
    }

    // H5Pclose(famFAPL_id);
    H5Pclose(fileAccessPropertiesList);
    return fileId;
}

hid_t OpenHdfFileReadOnly(std::string& fileName)
{
    hid_t fileAccessPropertiesList = CreateFAPL();

    hid_t fileId =
        H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, fileAccessPropertiesList);

    // H5Pclose(famFAPL_id);
    H5Pclose(fileAccessPropertiesList);
    return fileId;
}

hid_t OpenHdfFileReadOnlySerial(std::string& fileName)
{
    hid_t fileAccessPropertiesList = CreateSerialFAPL();

    hid_t fileId =
        H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, fileAccessPropertiesList);

    // H5Pclose(famFAPL_id);
    H5Pclose(fileAccessPropertiesList);
    return fileId;
}

hid_t OpenHdfFileSerial(std::string& fileName)
{
    MPI_Info info = MPI_INFO_NULL;
    // hid_t fileAccessPropertiesList = H5Pcreate(H5P_FILE_ACCESS);
    // H5Pset_fapl_mpio(fileAccessPropertiesList, MPI_COMM_WORLD, info);
    hid_t fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // H5Pclose(fileAccessPropertiesList);
    return fileId;
}

int RetrieveNumberOfLevelsFromFile(hid_t file_id)
{
    herr_t status;
    string numLevelsDataName = NUM_LEVELS_DATASET;
    hid_t numLevelsDataId =
        H5Dopen2(file_id, numLevelsDataName.c_str(), H5P_DEFAULT);
    int numLevels = -1;

    status = H5Dread(numLevelsDataId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &numLevels);
    if (numLevels == -1)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        cout << "Process rank " << rank
             << "did not import number of levels from file";
        MPI_Abort(MPI_COMM_WORLD, 1);
        MPI_Finalize();
        exit(1);
    }
    return numLevels;
}

int RetrieveNumberOfLevelsFromFileSerial(hid_t file_id)
{
    herr_t status;
    string numLevelsDataName = NUM_LEVELS_DATASET;
    hid_t numLevelsDataId =
        H5Dopen2(file_id, numLevelsDataName.c_str(), H5P_DEFAULT);
    int numLevels = -1;

    status = H5Dread(numLevelsDataId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, &numLevels);
    if (numLevels == -1 || status != 0)
    {
        cout << "did not import number of levels from file";

        exit(1);
    }
    return numLevels;
}

hid_t CreateFAPL()
{
    MPI_Info info;
    MPI_Info_create(&info);

    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#ifdef WIN32

    /* Disables ROMIO's data-sieving */
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");

    /* Enable ROMIO's collective buffering */
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");
#else

    /* Disables ROMIO's data-sieving */
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");

    /* Enable ROMIO's collective buffering */
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");

#endif
#ifdef CLUSTER

    long long stripeSize = 1024LL * 1024LL * 1024LL;
    string stripeString = std::to_string(stripeSize);
    MPI_Info_set(info, "cb_buffer_size", stripeString.c_str());
    MPI_Info_set(info, "striping_unit", stripeString.c_str());
    long long size = mpiSize;
    string numStripes = std::to_string(size);
    if (mpiSize < 24)
    {
        long long size = mpiSize;
        string numStripes = std::to_string(size);
        MPI_Info_set(info, "cb_nodes", numStripes.c_str());
        MPI_Info_set(info, "striping_factor", numStripes.c_str());
    }
    else
    {
        string numStripes = std::to_string(16LL);
        MPI_Info_set(info, "cb_nodes", numStripes.c_str());
        MPI_Info_set(info, "striping_factor", "24");
    }
#endif

    hid_t memFaplId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(memFaplId, MPI_COMM_WORLD, info);
    // hid_t famFaplId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(memFaplId, H5F_CLOSE_STRONG);
    H5Pset_libver_bounds(memFaplId, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    // H5Pset_meta_block_size(memFaplId, 100 * 1024 * 1024);
    H5Pset_meta_block_size(memFaplId, 32 * 1024 * 1024);
    // H5Pset_sieve_buf_size(memFaplId, 1024 * 1024 * 32);
    // H5Pset_cache(memFaplId, 2, 6421, 32 * 1024 * 1024, 0.8);
    H5Pset_small_data_block_size(memFaplId, 32 * 1024 * 1024);
#ifdef WIN32
    H5Pset_alignment(memFaplId, 16384, 4096);
#endif    // WIN32

#ifdef CLUSTER
    H5Pset_alignment(memFaplId, 8 * 16384, 4096);
    // H5Pset_alignment(memFaplId, 16384, 4096);
    // H5Pset_alignment(memFaplId, 1024*1024, 1024);
    // H5Pset_alignment(memFaplId, 1024*1024, 1024);
    // H5Pset_page_buffer_size(memFaplId, 1024 * 1024 * 1024, 0, 0);
    // H5Pset_file_space_strategy(memFaplId, H5F_FSPACE_STRATEGY_AGGR, true,
    // 1024);

#endif    // WIN32

    H5AC_cache_config_t mdc_config;
    // mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    // H5Pget_mdc_config(memFaplId, &mdc_config);

    // mdc_config.evictions_enabled = false;
    // mdc_config.incr_mode = H5C_incr__off;
    // mdc_config.decr_mode = H5C_decr__off;
    // mdc_config.flash_incr_mode = H5C_flash_incr__off;
    ////mdc_config.metadata_write_strategy =
    ///H5AC_METADATA_WRITE_STRATEGY__DISTRIBUTED;
    // H5Pset_mdc_config(memFaplId, &mdc_config);
    // H5Pset_all_coll_metadata_ops(memFaplId, true);

    return memFaplId;
}

hid_t CreateSerialFAPL()
{
    hid_t memFaplId = H5Pcreate(H5P_FILE_ACCESS);
    // hid_t famFaplId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(memFaplId, H5F_CLOSE_STRONG);
    H5Pset_libver_bounds(memFaplId, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    // H5Pset_meta_block_size(memFaplId, 100 * 1024 * 1024);
    H5Pset_meta_block_size(memFaplId, 32 * 1024 * 1024);
    // H5Pset_sieve_buf_size(memFaplId, 1024 * 1024 * 32);
    // H5Pset_cache(memFaplId, 2, 6421, 32 * 1024 * 1024, 0.8);
    H5Pset_small_data_block_size(memFaplId, 32 * 1024 * 1024);
#ifdef WIN32
    H5Pset_alignment(memFaplId, 16384, 4096);
#endif    // WIN32

#ifdef CLUSTER
    H5Pset_alignment(memFaplId, 8 * 16384, 4096);
    // H5Pset_alignment(memFaplId, 16384, 4096);
    // H5Pset_alignment(memFaplId, 1024*1024, 1024);
    // H5Pset_alignment(memFaplId, 1024*1024, 1024);
    // H5Pset_page_buffer_size(memFaplId, 1024 * 1024 * 1024, 0, 0);
    // H5Pset_file_space_strategy(memFaplId, H5F_FSPACE_STRATEGY_AGGR, true,
    // 1024);

#endif    // WIN32

    H5AC_cache_config_t mdc_config;
    // mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    // H5Pget_mdc_config(memFaplId, &mdc_config);

    // mdc_config.evictions_enabled = false;
    // mdc_config.incr_mode = H5C_incr__off;
    // mdc_config.decr_mode = H5C_decr__off;
    // mdc_config.flash_incr_mode = H5C_flash_incr__off;
    ////mdc_config.metadata_write_strategy =
    ///H5AC_METADATA_WRITE_STRATEGY__DISTRIBUTED;
    // H5Pset_mdc_config(memFaplId, &mdc_config);
    // H5Pset_all_coll_metadata_ops(memFaplId, true);

    return memFaplId;
}

hid_t CreateHdfFile(std::string filename)
{
    hid_t faplId = CreateFAPL();

    hid_t fcplId = H5Pcreate(H5P_FILE_CREATE);
    H5Pset_sym_k(fcplId, 16, 8);
    // H5Pset_file_space_page_size(fcplId, 4096);
    // H5Pset_istore_k(fcplId, 32);

    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, fcplId, faplId);

    if (file_id < 0)
    {
        cout << "File was not created" << file_id << endl;
        MPI_Abort(MPI_COMM_WORLD, file_id);
    }

    H5Pclose(faplId);
    H5Pclose(fcplId);

    // abort();

    return file_id;
}

hid_t CreateHdfFileSerial(std::string filename)
{
    hid_t faplId = CreateSerialFAPL();

    hid_t fcplId = H5Pcreate(H5P_FILE_CREATE);
    H5Pset_sym_k(fcplId, 16, 8);
    // H5Pset_file_space_page_size(fcplId, 4096);
    // H5Pset_istore_k(fcplId, 32);

    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, fcplId, faplId);

    if (file_id < 0)
    {
        cout << "File was not created" << file_id << endl;
        MPI_Abort(MPI_COMM_WORLD, file_id);
    }

    H5Pclose(faplId);
    H5Pclose(fcplId);

    // abort();

    return file_id;
}

void ExportConnectedVesselData(hid_t fileId, long long*& connectedVesselData)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t datasetId =
        OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    hsize_t numElemetsPerProc = dims[0] / mpiSize;
    hsize_t localNumElements = numElemetsPerProc;
    if (mpiRank == mpiSize - 1)
    {
        localNumElements += dims[0] % mpiSize;
    }

    // cout << "Rank " << mpiRank << " Num Elements: " << localNumElements <<
    // endl;

    hsize_t start[2], block[2], count[2], stride[2];

    stride[0] = 1;
    stride[1] = 1;
    block[0] = 1;
    block[1] = 1;
    start[0] = numElemetsPerProc * mpiRank;
    start[1] = 0;
    count[0] = localNumElements;
    count[1] = 1;
    // cout << " MPI Rank: " << mpiRank << " Start: " << start[0] << " Count: " <<
    // count[0] << endl;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);

    hid_t dataTransferPropertiesList = H5Pcreate(H5P_DATASET_XFER);
    if (mpiSize > 1)
        H5Pset_dxpl_mpio(dataTransferPropertiesList, H5FD_MPIO_COLLECTIVE);

    hid_t memSpaceId = H5Screate_simple(1, &localNumElements, NULL);

    H5Dwrite(datasetId, H5T_NATIVE_LLONG, memSpaceId, dataspace,
        dataTransferPropertiesList, connectedVesselData);

    H5Pclose(dataTransferPropertiesList);
    H5Dclose(datasetId);
    H5Sclose(dataspace);
    H5Sclose(memSpaceId);
}

void ExportConnectedVesselData(
    hid_t fileId, std::vector<long long>& connectedVesselData)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t datasetId =
        OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);

    hsize_t numElemetsPerProc = GetNumRowsInDataset(datasetId) / mpiSize;
    hsize_t localNumElements = connectedVesselData.size();

    hsize_t localStartIndex = numElemetsPerProc * mpiRank;

    long long remainingNumberOfVesselsToExport =
        localNumElements % (maxNumVesselsPerWrite * 4);
    int numLoopsNeeded = localNumElements / (maxNumVesselsPerWrite * 4);
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(
        remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    if (remainingNumberOfVesselsToExport > 0)
    {
        CollectiveHdfBlockExport(datasetId, localStartIndex,
            remainingNumberOfVesselsToExport, connectedVesselData, 0);
    }
    else
    {
        CollectiveHdfNullImport<long long>(datasetId);
    }
    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
            hsize_t startIndex = localStartIndex +
                remainingNumberOfVesselsToExport +
                (4 * maxNumVesselsPerWrite * i);
            hsize_t dataOffset = remainingNumberOfVesselsToExport +
                (4 * maxNumVesselsPerWrite * i);

            CollectiveHdfBlockExport(datasetId, startIndex,
                maxNumVesselsPerWrite * 4, connectedVesselData, dataOffset);
        }
        else
        {
            CollectiveHdfNullExport<long long>(datasetId);
        }
    }

    H5Dclose(datasetId);
}

long long ImportConnectedVesselData(
    hid_t fileId, long long*& connectedVesselData)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t datasetId =
        OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    hsize_t numElemetsPerProc = dims[0] / mpiSize;
    hsize_t localNumElements = numElemetsPerProc;
    if (mpiRank == mpiSize - 1)
    {
        localNumElements += dims[0] % mpiSize;
    }

    // cout << "Rank " << mpiRank << " Num Elements: " << localNumElements <<
    // endl;

    hsize_t start[2], block[2], count[2], stride[2];

    stride[0] = 1;
    stride[1] = 1;
    block[0] = 1;
    block[1] = 1;
    start[0] = numElemetsPerProc * mpiRank;
    start[1] = 0;
    count[0] = localNumElements;
    count[1] = 1;
    // cout << " MPI Rank: " << mpiRank << " Start: " << start[0] << " Count: " <<
    // count[0] << endl;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, block);

    connectedVesselData = new long long[localNumElements];

    hid_t dataTransferPropertiesList = H5Pcreate(H5P_DATASET_XFER);
    if (mpiSize > 1)
        H5Pset_dxpl_mpio(dataTransferPropertiesList, H5FD_MPIO_COLLECTIVE);

    hid_t memSpaceId = H5Screate_simple(1, &localNumElements, NULL);

    H5Dread(datasetId, H5T_NATIVE_LLONG, memSpaceId, dataspace,
        dataTransferPropertiesList, connectedVesselData);

    H5Pclose(dataTransferPropertiesList);
    H5Dclose(datasetId);
    H5Sclose(dataspace);
    H5Sclose(memSpaceId);

    return localNumElements;
}

std::vector<long long> ImportConnectedVesselData(hid_t fileId)
{
    int mpiSize, mpiRank;

    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    hid_t datasetId =
        OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    hsize_t numElemetsPerProc = dims[0] / mpiSize;
    hsize_t localNumElements = numElemetsPerProc;
    if (mpiRank == mpiSize - 1)
    {
        localNumElements += dims[0] % mpiSize;
    }
    hsize_t localStartIndex = numElemetsPerProc * mpiRank;

    long long remainingNumberOfVesselsToExport =
        localNumElements % (maxNumVesselsPerWrite * 4);
    int numLoopsNeeded = localNumElements / (maxNumVesselsPerWrite * 4);
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(
        remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    std::vector<long long> connectedVesselData;
    connectedVesselData.resize(localNumElements);
    if (remainingNumberOfVesselsToExport > 0)
    {
        CollectiveHdfAppendingImport(datasetId, localStartIndex,
            remainingNumberOfVesselsToExport, connectedVesselData, 0);
    }
    else
    {
        CollectiveHdfNullImport<long long>(datasetId);
    }
    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
            hsize_t startIndex = localStartIndex +
                remainingNumberOfVesselsToExport +
                (4 * maxNumVesselsPerWrite * i);
            hsize_t dataOffset = remainingNumberOfVesselsToExport +
                (4 * maxNumVesselsPerWrite * i);
            CollectiveHdfAppendingImport(datasetId, startIndex,
                maxNumVesselsPerWrite * 4, connectedVesselData, dataOffset);
        }
        else
        {
            CollectiveHdfNullImport<long long>(datasetId);
        }
    }

    H5Dclose(datasetId);
    return std::move(connectedVesselData);
}

int GetNumPartitionsInFile(hid_t fileId)
{
    int fileNumPartitions = 0;

    // hid_t ftplId = CreateCollectiveDataTransferPropertiesList();

    hid_t attributeId = H5Aopen(fileId, "Num_Partitions", H5P_DEFAULT);

    H5Aread(attributeId, H5T_NATIVE_INT, &fileNumPartitions);

    H5Aclose(attributeId);
    // H5Pclose(ftplId);
    return fileNumPartitions;
}

long long GetNumRowsInDataset(hid_t datasetId)
{
    hsize_t dims[2] = {0};
    hid_t dataspaceId = H5Dget_space(datasetId);
    H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
    H5Sclose(dataspaceId);
    return dims[0];
}

void ImportGeomDataBlock(hid_t datasetId, long long startVesselIndex,
    long long numVesselsToImport, std::vector<double>& data)
{
    CollectiveHdfBlockImport<double, 2, 7>(
        datasetId, startVesselIndex, numVesselsToImport, data);
    // std::vector<double> test(numVesselsToImport * 7 + 10);
    // hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    // hid_t dataspaceId = H5Dget_space(datasetId);
    // hsize_t start[2] = { 0 }, count[2] = { 0 }, strideAndBlocks[2] = { 1,1 };
    // count[0] = numVesselsToImport;
    // count[1] = 7;
    // start[0] = startVesselIndex;
    // start[1] = 0;
    // hid_t memspaceId = H5Screate_simple(2, count, NULL);
    // H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
    // count, strideAndBlocks);

    ////H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId,
    ///H5P_DEFAULT, data.data());
    // H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, dtplId,
    // test.data());

    // H5Sclose(dataspaceId);
    // H5Sclose(memspaceId);
    // H5Pclose(dtplId);

    // if (std::equal(test.begin(), test.end(), data.begin()))
    //{
    //	cout << "Geom Vectors are same" <<endl;
    //}
    // else
    //	cout << "Geom Vectors not the same" <<endl;
}

void ImportNodeDataBlock(hid_t datasetId, long long startVesselIndex,
    long long numVesselsToImport, std::vector<long long>& data)
{
    CollectiveHdfBlockImport<long long, 2, 2>(
        datasetId, startVesselIndex, numVesselsToImport, data);
    // std::vector<long long> test(numVesselsToImport * 2 + 10);
    // hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    // hid_t dataspaceId = H5Dget_space(datasetId);
    // hsize_t start[2] = { 0 }, count[2] = { 0 }, strideAndBlocks[2] = { 1,1 };
    // count[0] = numVesselsToImport;
    // count[1] = 2;
    // start[0] = startVesselIndex;
    // start[1] = 0;
    // hid_t memspaceId = H5Screate_simple(2, count, NULL);

    // H5Sselect_hyperslab(dataspaceId, H5S_SELECT_SET, start, strideAndBlocks,
    // count, strideAndBlocks);
    ////H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId,
    ///H5P_DEFAULT, data.data());
    // H5Dread(datasetId, H5T_NATIVE_LLONG, memspaceId, dataspaceId, dtplId,
    // test.data());

    // H5Sclose(dataspaceId);
    // H5Pclose(dtplId);
    // H5Sclose(memspaceId);
    // if (std::equal(test.begin(), test.end(), data.begin()))
    //{
    //	cout << "node Vectors are same" << endl;
    //}
    // else
    //	cout << "node Vectors not the same" << endl;
}

void ImportGeometryAndNodeDataBlock(hid_t geomDatasetId, hid_t nodeDatasetId,
    long long startVesselIndex, long long numVesselsToImport,
    std::vector<double>& geomData, std::vector<long long>& nodeData)
{
    ImportGeomDataBlock(
        geomDatasetId, startVesselIndex, numVesselsToImport, geomData);
    ImportNodeDataBlock(
        nodeDatasetId, startVesselIndex, numVesselsToImport, nodeData);
}

void NullImportGeometry(hid_t datasetId)
{
    double* nullDataPointer = nullptr;
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = 0;
    count[1] = 7;
    start[0] = 0;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(2, count, NULL);
    hid_t dataspaceId = H5Screate_simple(2, count, NULL);

    // H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, H5P_DEFAULT,
    // data.data());
    H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, dtplId,
        nullDataPointer);

    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

void NullImportNode(hid_t datasetId)
{
    double* nullDataPointer = nullptr;
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();

    hsize_t start[2] = {0}, count[2] = {0}, strideAndBlocks[2] = {1, 1};
    count[0] = 0;
    count[1] = 2;
    start[0] = 0;
    start[1] = 0;
    hid_t memspaceId = H5Screate_simple(2, count, NULL);
    hid_t dataspaceId = H5Screate_simple(2, count, NULL);

    // H5Dread(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, H5P_DEFAULT,
    // data.data());
    H5Dread(datasetId, H5T_NATIVE_LLONG, memspaceId, dataspaceId, dtplId,
        nullDataPointer);

    H5Sclose(dataspaceId);
    H5Sclose(memspaceId);
    H5Pclose(dtplId);
}

void PerformNullImportOfGeometryAndNodes(
    hid_t geomDatasetId, hid_t nodeDatasetId)
{
    CollectiveHdfNullImport<double, 2, 7>(geomDatasetId);
    CollectiveHdfNullImport<long long, 2, 2>(nodeDatasetId);
}

void DetermineIOLoopCount(long long& remainingNumberOfVesselsToExport,
    int& numLoopsNeeded, int& globalNumLoopsNeeded)
{
    int mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    if (mpiSize > 1)
    {
        bool localNumRemainingVessels = remainingNumberOfVesselsToExport;
        bool anyProcessHasRemainingVessels = false;
        MPI_Allreduce(&localNumRemainingVessels, &anyProcessHasRemainingVessels,
            1, MPI_CHAR, MPI_LOR, MPI_COMM_WORLD);
        if (anyProcessHasRemainingVessels && !localNumRemainingVessels &&
            numLoopsNeeded > 0)
        {
            remainingNumberOfVesselsToExport = maxNumVesselsPerWrite;
            numLoopsNeeded -= 1;
        }

        MPI_Allreduce(&numLoopsNeeded, &globalNumLoopsNeeded, 1, MPI_INT,
            MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
        globalNumLoopsNeeded = numLoopsNeeded;
    }
}

std::vector<long long> ImportNumberOfVesselsPerPartitionDataOld(
    hid_t fileId, int numPartitions)
{
    hid_t datasetId =
        OpenHdfDataset(fileId, "/", NUM_VESSELS_PER_PARTITION_DATASET);
    hid_t dataspace = H5Dget_space(datasetId);
    hsize_t dims[2] = {0};
    // hsize_t numDims = H5Sget_simple_extent_ndims(dataspace);
    // cout << "Num Dims: " << numDims;
    H5Sget_simple_extent_dims(dataspace, dims, NULL);

    // long long*  numVesselsInPartitions = new long long[numPartitions];
    std::vector<long long> numVesselsPerPartition(numPartitions);
    hid_t dtplId = CreateCollectiveDataTransferPropertiesList();
    H5Dread(datasetId, H5T_NATIVE_LLONG, H5S_ALL, dataspace, dtplId,
        numVesselsPerPartition.data());
    H5Sclose(dataspace);
    H5Dclose(datasetId);
    H5Pclose(dtplId);

    return numVesselsPerPartition;
}

std::vector<long long> ImportNumberOfVesselsPerPartitionData(
    hid_t fileId, int numPartitions)
{
    std::vector<long long> numVesselsPerPartition(numPartitions);
    hid_t datasetId =
        OpenHdfDataset(fileId, "/", NUM_VESSELS_PER_PARTITION_DATASET);
    CollectiveHdfBlockImport(
        datasetId, 0, numPartitions, numVesselsPerPartition);
    H5Dclose(datasetId);
    return numVesselsPerPartition;
}

long long GetTotalNumberOfVessels(hid_t fileId)
{
    hid_t geomDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, GEOM_DATASET);
    hsize_t dims[2] = {0};
    hid_t dataspaceId = H5Dget_space(geomDatasetId);
    H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
    H5Sclose(dataspaceId);
    H5Dclose(geomDatasetId);
    return dims[0];
}