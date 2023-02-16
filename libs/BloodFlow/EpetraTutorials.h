#pragma once

// Trilinos Epetra headers
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include "VesselGenerator/network.h"
#include <string>

void exampleRoutine (const Epetra_Comm& comm, std::ostream& out);

double powerMethod(const Epetra_Operator&A, const int niters, const double tolerance);

void EPetraCrsMatrix();

void AssembleMatrixSequential(Network& network, std::string filename, int numLevels);

void MLAztecOO();

bool CheckFlowsAreCorrect(Network& network, Epetra_Vector& flows);