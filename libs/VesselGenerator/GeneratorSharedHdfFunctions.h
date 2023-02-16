#pragma once

#include <vector>

#include <hdf5.h>

#include "vessel.h"

typedef std::vector<Vessel> VesselVector;

/**
 * \breif Function to serialize vessel data for root vessels in parallel algorithm
 * Size of geomArray must be 7 times the number of vessels
 * Size of nodeArray must be 7 times the number of vessels
 * @param[in] numVesselsToWrite
 * @param[in] arteries
 * @param[in] veins
 * @param[in,out] geomArray
 * @param[in,out] nodeArray
 */

void ConvertRootVesselVectorsToDataTables(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, double* geomArray,
    long long* nodeArray, double* conductanceArray);
/**
 * Function to serialize vessel objects
 *
 * Size of geomArray must be 7 times the number of vessels
 * Size of nodeArray must be 7 times the number of vessels
 * @param[in] numVesselsToWrite
 * @param[in] arteries
 * @param[in] veins
 * @param[in,out] geomArray
 * @param[in,out] nodeArray
 */
void ConvertVesselVectorsToDataTables(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, double* geomArray,
    long long* nodeArray);
/**
 * Function that serializes the start pt, end pt, and radius of a vessel object
 *
 * Size of geomArray must be 7 times the number of vessels
 * @param[in] numVesselsToWrite
 * @param[in] vessels
 * @param[out] geomArray
 */
void ConvertVesselVectorsToGeometryArray(long long numVesselsToWrite,
    VesselVector& vessels, double* geomArray);

/**
 * Function to serialize node identifiers in vessel objects
 * Size of nodeArray must be 2 times the number of vessels
 * @param numVesselsToWrite
 * @param arteries
 * @param veins
 * @param nodeArray
 */
void ConvertVesselVectorsToNodeArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* nodeArray);

/**
 * Function to serialize conductances in vessel objects
 * Size of conductanceArray must be the number of vessels
 * @param numVesselsToWrite
 * @param arteries
 * @param veins
 * @param conductanceArray
 */
void ConvertVesselVectorsToConductanceArray(long long numVesselsToWrite,
    VesselVector& arteries, VesselVector& veins, double* conductanceArray);

/**
 * Creates all groups in HDF5 file
 * @param file_id
 */
void CreateHdfGroups(hid_t file_id);

/**
 * Create HDF5 dataset and populate it with number of generations in network
 * @param file_id
 * @param levels
 */
void CreateAndExportNumberOfLevels(hid_t file_id, int levels);

/**
 * Create a dataset with the bounding box in the HDF5 file.
 * @param file_id
 * @param boundingBox
 */
void CreateAndExportBoundingBox(hid_t file_id, std::vector<double>& boundingBox);

/**
 * Create the HDF5 dataset for the geometry data
 * @param file_id
 * @param levels
 * @param numVessels
 * @return dataset handle
 */
hid_t CreateGeometryHdfDataset(hid_t file_id, int levels, long long numVessels);

/**
 * Create the HDF5 dataset for the node data
 * @param file_id
 * @param levels
 * @param numVessels
 * @return dataset handle
 */
hid_t CreateNodeHdfDataset(hid_t file_id, int levels, long long numVessels);

/**
 * Create the HDF5 dataset for the healthy conductance data
 * @param file_id
 * @param levels
 * @param numVessels
 * @return dataset handle
 */
hid_t CreateHealthyConductanceHdfDataset(hid_t file_id, int levels, long long numVessels);

/**
 * Function to serialize the number of connected vessels for each vessel in network
 *
 * numConnectedVesselsArray must be equal in length to the number of vessels
 * @param arteries
 * @param veins
 * @param numConnectedVesselsArray
 * @param numVessels
 * @return
 */
long long PopulateRootNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins, short* numConnectedVesselsArray,
    long long numVessels);

/**
 * Function to serialize the number of vessel connections for the root network in the parallel algorithm
 *
 * numConnectedVesselsArray must be equal in length to the number of vessels
 * @param arteries
 * @param veins
 * @param numConnectedVesselsArray
 * @param numVessels
 * @return total number of connected vessels
 */
long long PopulateNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins, short* numConnectedVesselsArray,
    long long numVessels);

/**
 * Serialize the indexes of connected vessels for the vessels into an array
 *
 * connectedVesselArray must be initialized to the total number of connected vessels
 * @param numVesselsToWrite
 * @param arteries
 * @param veins
 * @param connectedVesselArray
 */
void PopulateRootConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* connectedVesselArray);

/**
 * Serialize the indexes of connected vessels for the vessels into an array
 *
 * connectedVesselArray must be initialized to the total number of connected vessels
 * @param numVesselsToWrite
 * @param arteries
 * @param veins
 * @param connectedVesselArray
 */
void PopulateConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins, long long* connectedVesselArray);

/**
 * Creates attributes to describe the number of partitions in the network
 * @param fileId
 */
void CreatePartitionAttributes(hid_t fileId);

/**
 * Creates the dataset for the number of connected vessels
 * @param file_id
 * @param levels
 * @param numVessels
 * @return
 */
hid_t createNumConnectedVesselsDataset(
    hid_t file_id, int levels, long long numVessels);

/**
 * Creates dataset for the connected vessel indexes
 * @param file_id
 * @param levels
 * @param totalNumConnectedVessels
 * @return
 */
hid_t createConnectedVesselsDataset(hid_t file_id, int levels,
    std::vector<unsigned long>::size_type totalNumConnectedVessels);

/**
 * Function to export the blood flow boundary conditions for the network
 *
 * This should probably be used as command-line arguments in either this or the blood flow algorithm.
 * @param file_id
 * @param levels
 */
void CreateAndExportBoundaryConditions(hid_t file_id, int levels);

/**
 * Vector based function for serialize the number of connected vessels
 * @param arteries
 * @param veins
 * @param numConnectedVesselsArray
 * @param numVessels
 * @return
 */
long long PopulateNumConnectedVesselsArray(
    VesselVector& arteries,
    VesselVector& veins,
    std::vector<short>& numConnectedVesselsArray, long long numVessels);

/**
 * Vector based version of algorithm that serializes the indexes of connected vessels.
 * @param numVesselsToWrite
 * @param arteries
 * @param veins
 * @param connectedVesselArray
 */
void PopulateConnectedVesselArray(long long numVesselsToWrite,
    VesselVector& arteries,
    VesselVector& veins,
    std::vector<long long>& connectedVesselArray);