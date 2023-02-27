#pragma once

#include <string>
#include <vector>

#include <hdf5.h>

#include "PreprocessorVessel.h"

/**
 * Function that exports the sorted vessel data back to the HDF5 file
 * @param fileId
 * @param vessels
 * @param numParts
 */
void ExportVesselDataSorted(hid_t fileId, std::vector<PreprocessorVessel>& vessels, int numParts);

/**
 * Function that calculates where the vessels on each compute node must be written to the HDF5 file
 * @param vessels
 * @param numParts
 * @param localNumVesselsPerPartition
 * @param outputOffsets
 */
void CalculateVesselOffsets(std::vector<PreprocessorVessel>& vessels, int numParts, long long*& localNumVesselsPerPartition, hsize_t*& outputOffsets);

/**
 * Function that calculates where the connected vessel data on each compute node must be written to the HDF5 file
 * @param vessels
 * @param numParts
 * @param outputOffsets
 */
// void CalculateConnectedVesselOffsets(std::vector<PreprocessorVessel>& vessels, int numParts, hsize_t*& outputOffsets);

/**
 * Function that calculates where the connected vessel data on each compute node must be written to the HDF5 file
 * @param vessels
 * @param numParts
 * @param outputOffsets
 */
void CalculateConnectedVesselOffsets(std::vector<PreprocessorVessel>& vessels, int numParts, std::vector<hsize_t>& outputOffsets);