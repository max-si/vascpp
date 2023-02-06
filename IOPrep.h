#pragma once

#include <vector>

#include "vessel.h"

typedef std::vector<Vessel> VesselVector;

// *******
void ExportRootVesselsToHdf5_HF(std::string filename, int totalLevels, int rootLevels,
    VesselVector& arteries,
    VesselVector& veins);

// *******
void ExportSubVesselsToHdf5_HF(std::string filename, int totalLevels, int subLevels,
    VesselVector& arteries, VesselVector& veins);

//*******
void ExportSequential_HF(std::string filename, int totalLevels, VesselVector& vessels);

//*******
void prepFile_HF(std::string filename, int numLevels, std::vector<double>& bbox, double output_bc);