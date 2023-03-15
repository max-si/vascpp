#pragma once

// Trilinos Epetra headers
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include "VesselGenerator/network.h"
#include <string>

void AssembleMatrixSequential(Network& network, std::string filename, int numLevels);

void MLAztecOO();

bool CheckFlowsAreCorrect(Network& network, Epetra_Vector& flows);

void ImportNodesAndConductances(std::string filename, int start, int numRows, std::vector<long long>& nodeArray, std::vector<double>& conductanceArray);