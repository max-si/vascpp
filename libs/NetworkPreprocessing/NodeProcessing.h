#pragma once

#include <vector>

#include <hdf5.h>

/**
 * Function that specifies the edge junction finding algorithm
 * This function is a delf contaitned unit.
 * It assumes that the vessels in the file have been correctly grouped by the partition and the vessel indexes
 * have been updated to represent new vessel locations.
 *
 * @param fileId
 */
void NodeProcessing(hid_t fileId);