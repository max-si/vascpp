#pragma once

#include <hdf5.h>

#include "NetworkDescription.h"

/**
 *  Function that exports blood flow rates to the HDF5 file
 *
 * @param fileId HDF5 file handle to write to
 * @param network Network Description of vessel data to Write to
 */
void BloodFlowExporter(hid_t fileId, const NetworkDescription& network);