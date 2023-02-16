
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <string>

#include <mpi.h>
#include <hdf5.h>

#include "PreprocessorVessel.h"
#include "preImporter.h"
#include "VesselGenerator/coordinate.h"
#include "VesselGenerator/ImportExportCommon.h"
#include "core/Timsort.h"

typedef std::vector<PreprocessorVessel> VesselVector;
typedef unsigned long ulong;

void ConstructVesselsFromData(long long numVesselsToConstruct,
    std::vector<PreprocessorVessel>& vessels, std::vector<double>& geomData,
    std::vector<long long>& nodeData)
{
    for (long long i = 0; i < numVesselsToConstruct; i++)
    {
        //cout << i << endl;
        std::vector<double>::size_type geomIndex = i * 7;
        std::vector<double>::size_type nodeIndex = i * 2;
        Coordinate startPoint(geomData[geomIndex], geomData[geomIndex + 1], geomData[geomIndex + 2]);
        Coordinate endPoint(geomData[geomIndex + 3], geomData[geomIndex + 4], geomData[geomIndex + 5]);
        double radius = geomData[geomIndex + 6];
        long long node1 = nodeData[nodeIndex];
        long long node2 = nodeData[nodeIndex + 1];
        PreprocessorVessel currentVessel(startPoint, endPoint, radius, node1, node2);
        vessels.push_back(currentVessel);
    }
}

std::vector<PreprocessorVessel> ImportGeometryAndConstructBaseVessels(hid_t fileId)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    
    hid_t geomDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, GEOM_DATASET);
    hid_t nodeDatasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NODE_DATASET);

    hsize_t dims[2] = {0};
    hid_t dataspaceId = H5Dget_space(geomDatasetId);
    H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
    H5Sclose(dataspaceId);
    long long localNumVessels = dims[0] / mpiSize;

    long long vesselOffset = localNumVessels * mpiRank;

    if (mpiRank == mpiSize - 1) { localNumVessels += dims[0] % mpiSize; }

    // determine import block size
    long long remainingNumberOfVesselsToExport = localNumVessels % maxNumVesselsPerWrite;
    int numLoopsNeeded = localNumVessels / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    //createRequiredArrays
    std::vector<double> geomImportArray;
    geomImportArray.resize(maxNumVesselsPerWrite * 7);
    std::vector<double> geomBuildArray;
    geomBuildArray.resize(maxNumVesselsPerWrite * 7);

    std::vector<long long> nodeImportArray;
    nodeImportArray.resize(maxNumVesselsPerWrite * 2);
    std::vector<long long> nodeBuildArray;
    nodeBuildArray.resize(maxNumVesselsPerWrite * 2);
    std::vector<PreprocessorVessel> vessels;
    vessels.reserve(localNumVessels);

    long long numVesselsToConstruct = 0;
    if (remainingNumberOfVesselsToExport > 0)
    {
        ImportGeometryAndNodeDataBlock(geomDatasetId, nodeDatasetId,
            vesselOffset, remainingNumberOfVesselsToExport, geomBuildArray, nodeBuildArray);
        numVesselsToConstruct = remainingNumberOfVesselsToExport;
    }

    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
            long long vesselStartLocation = vesselOffset + remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;
//#pragma omp parallel sections num_threads(2)
            {
//#pragma omp section
                {
                    ImportGeometryAndNodeDataBlock(geomDatasetId, nodeDatasetId,
                        vesselStartLocation, maxNumVesselsPerWrite,
                        geomImportArray, nodeImportArray);
                }
//#pragma omp section
                {
                    ConstructVesselsFromData(numVesselsToConstruct, vessels, geomBuildArray, nodeBuildArray);
                }
            }

            std::swap(geomBuildArray, geomImportArray);
            std::swap(nodeBuildArray, nodeImportArray);
            numVesselsToConstruct = maxNumVesselsPerWrite;
        }
        else
        {
            PerformNullImportOfGeometryAndNodes(geomDatasetId, nodeDatasetId);
        }
    }
    ConstructVesselsFromData(numVesselsToConstruct, vessels, geomBuildArray, nodeBuildArray);
//#pragma omp parallel for default(shared)
    for (long long i = 0; i < localNumVessels; ++i)
    {
        vessels[i].setInitialArrayPosition(vesselOffset + i);
    }
    H5Dclose(geomDatasetId);
    H5Dclose(nodeDatasetId);
    return std::move(vessels);
}

std::vector<short> ImportNumConnectedVessels(hid_t fileId, long long localVesselCount, long long localStartIndex)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport =localVesselCount % maxNumVesselsPerWrite;
    int numLoopsNeeded = localVesselCount / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    std::vector<short> importData;
    importData.resize(localVesselCount);

    hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, NUM_CONNECTED_VESSELS_DATASET);

    if (remainingNumberOfVesselsToExport > 0)
    {
        CollectiveHdfAppendingImport(datasetId, localStartIndex, remainingNumberOfVesselsToExport, importData, 0);
    }

    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
            long long vesselStartLocation = localStartIndex +
                remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite;

            CollectiveHdfAppendingImport(datasetId, vesselStartLocation,
                maxNumVesselsPerWrite, importData,
                remainingNumberOfVesselsToExport + i * maxNumVesselsPerWrite);
        }
        else
        {
            CollectiveHdfNullImport<short>(datasetId);
        }
    }

    H5Dclose(datasetId);
    return importData;
}

void AddConnectedVesselsToVessels(std::vector<PreprocessorVessel>& vessels, std::vector<short>& numConnectedVesselsArray,
    std::vector<long long>& connectedVesselArray, long long startVesselIndex, long long numVesselsToConstruct)
{
    std::vector<long long>::iterator startItr = connectedVesselArray.begin();
    for (long long i = startVesselIndex;
         i < startVesselIndex + numVesselsToConstruct; ++i)
    {
        short numConnectedVessels = numConnectedVesselsArray[i];
        std::vector<long long>::iterator endItr = startItr + numConnectedVessels;
        vessels[i].assignConnectedVessels(startItr, endItr);
        startItr = endItr;
    }
}

void ImportConnectedVesselIndexesAndAddToVessels(hid_t fileId,
    std::vector<short>& numConnectedVesselsArray, long long localVesselCount,
    long long dataStartIndex, std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    long long remainingNumberOfVesselsToExport = localVesselCount % maxNumVesselsPerWrite;
    int numLoopsNeeded = localVesselCount / maxNumVesselsPerWrite;
    int globalNumLoopsNeeded = 0;

    DetermineIOLoopCount(remainingNumberOfVesselsToExport, numLoopsNeeded, globalNumLoopsNeeded);

    hid_t datasetId = OpenHdfDataset(fileId, GEOM_GROUP_NAME, CONNECTED_VESSELS_DATASET);

    //createRequiredArrays
    std::vector<long long> importArray;
    std::vector<long long> buildArray;

    long long numConnectedVesselsToImport;
    long long numVesselsToConstruct = 0;
    long long importStartIndex = dataStartIndex;
    long long constructionStartIndex = 0;
    if (remainingNumberOfVesselsToExport > 0)
    {
        numConnectedVesselsToImport = std::accumulate(
            numConnectedVesselsArray.begin(),
            numConnectedVesselsArray.begin() + remainingNumberOfVesselsToExport,
            0LL);
        CollectiveHdfBlockImport(datasetId, importStartIndex,
            numConnectedVesselsToImport, buildArray);
        importStartIndex += numConnectedVesselsToImport;
        numVesselsToConstruct = remainingNumberOfVesselsToExport;
    }
    for (int i = 0; i < globalNumLoopsNeeded; i++)
    {
        if (i < numLoopsNeeded)
        {
//#pragma omp parallel sections num_threads(2)
            {
//#pragma omp section
                {
                    numConnectedVesselsToImport =
                        std::accumulate(numConnectedVesselsArray.begin() +
                                remainingNumberOfVesselsToExport +
                                (i) *maxNumVesselsPerWrite,
                            numConnectedVesselsArray.begin() +
                                remainingNumberOfVesselsToExport +
                                (i + 1) * maxNumVesselsPerWrite,
                            0LL);
                    CollectiveHdfBlockImport(datasetId, importStartIndex,
                        numConnectedVesselsToImport, importArray);
                    importStartIndex += numConnectedVesselsToImport;
                }
//#pragma omp section
                {
                    AddConnectedVesselsToVessels(vessels,
                        numConnectedVesselsArray, buildArray,
                        constructionStartIndex, numVesselsToConstruct);
                    constructionStartIndex += numVesselsToConstruct;
                }
            }

            std::swap(buildArray, importArray);
            numVesselsToConstruct = maxNumVesselsPerWrite;
        }
        else
        {
            CollectiveHdfNullImport<long long>(datasetId);
        }
    }

    AddConnectedVesselsToVessels(vessels, numConnectedVesselsArray, buildArray, constructionStartIndex, numVesselsToConstruct);
    H5Dclose(datasetId);
}

void ImportConnectedVesselData(hid_t fileId, std::vector<PreprocessorVessel>& vessels)
{
    int mpiSize, mpiRank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    long long totalNumVessels = GetTotalNumberOfVessels(fileId);
    long long localNumVessels = totalNumVessels / mpiSize;

    long long vesselOffset = localNumVessels * mpiRank;
    if (mpiRank == mpiSize - 1) { localNumVessels += totalNumVessels % mpiSize; }

    //import Number Connected Vessels
    std::vector<short> numConnectedVessels = ImportNumConnectedVessels(fileId, localNumVessels, vesselOffset);

    //Determine Import Offsets and Numbers
    long long localNumConnectedVessels = std::accumulate(numConnectedVessels.begin(), numConnectedVessels.end(), 0LL);
    long long startPosition = 0;
    MPI_Exscan(&localNumConnectedVessels, &startPosition, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpiRank == 0) { startPosition = 0; }

    //Import Connected Vessel Information
    ImportConnectedVesselIndexesAndAddToVessels(fileId, numConnectedVessels, localNumVessels, startPosition, vessels);
}

std::vector<PreprocessorVessel> PreprocessorVesselImporter(hid_t fileId) {
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	std::vector<PreprocessorVessel> vessels = ImportGeometryAndConstructBaseVessels(fileId);

	ImportConnectedVesselData(fileId, vessels);

	return std::move(vessels);
}