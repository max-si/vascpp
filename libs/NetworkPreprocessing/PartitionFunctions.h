#pragma once

#include <string>
#include <vector>

#include "AdjacencyMatrix.h"
#include "PreprocessorVessel.h"

// temp function to act as "main"
void networkPreproc(std::string filename);

// Function to create adjacency matrix for vessels
// Performs a communication on MPI_COMM_WORLD to collect network offsets based on the adjacency matrix
void GenerateAdjacencyMatrix(std::vector<PreprocessorVessel>& vessels, int numParts, AdjacencyMatrix& AdjacencyMatrix);

// Function that calls the ParMETIS algorithm on the adjacency matrix
void PartitionMatrix(std::vector<PreprocessorVessel>& vessels, AdjacencyMatrix& AdjacencyMatrix);

// In-place sorting of vessel objects to group local vessels by partition
// This facilitates renumbering and exporting
void SortVesselsByPartition(std::vector<PreprocessorVessel>& vessels);