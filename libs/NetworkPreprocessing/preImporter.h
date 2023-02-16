#pragma once

#include <string>
#include <vector>

#include <hdf5.h>

#include "PreprocessorVessel.h"

// parallel importer
//int PreprocessorVesselImporter(hid_t fileId, std::vector<PreprocessorVessel>& vessels);

std::vector<PreprocessorVessel> PreprocessorVesselImporter(hid_t fileId);