#pragma once

#include <vector>
#include <string>

#include <hdf5.h>

#include "VesselGenerator/vessel.h"
#include "Triplet.h"

typedef std::vector<Triplet> TripletVector;

void CreatePressureTripletList(std::string filename, int numLevels, TripletVector& TripletList);

void ImportNodesAndConductances(std::string filename, int start, int numRows, std::vector<long long>& nodeArray, std::vector<double>& conductanceArray);