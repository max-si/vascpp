#pragma once

#include <fstream>
#include <iostream>
#include <vector>

#include <hdf5.h>

#include "vessel.h"

typedef std::vector<Vessel> VesselVector;

/**
 * Function that exports the root tree to a new HDF5 file.
 *
 * This function creates the file, the datastructures, and exports the root tree.
 * @param filename
 * @param arteries
 * @param veins
 * @param levels
 * @param rootLevels
 * @param boundingBox
 * @return file handle
 */
hid_t RootVesselExporter(std::string filename,
    VesselVector& arteries,
    VesselVector& veins, int levels, int rootLevels,
    std::vector<double>& boundingBox);
/**
 * fuction that exports the sub-networks to the HDF5 file.
 * performs the export using collective properties.
 * @param file_id
 * @param arteries
 * @param veins
 * @param levels
 * @param subLevels
 */
void SubVesselExporter(hid_t file_id,
    VesselVector& arteries,
    VesselVector& veins, int levels, int subLevels);

/**
 * Serial exporter for complete network
 * @param filename
 * @param vessels
 * @param levels
 * @param boundingBox
 */
void VesselNetworkExporter_HDF5(std::string filename, int levels, VesselVector& vessels, std::vector<double>& boundingBox);