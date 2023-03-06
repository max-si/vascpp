#pragma once

#include <vector>

#include "NetworkDescription.h"

/**
 * Function to calculate healthy pressures in the network
 * @param network
 * @param pressureSolutionVector
 */
void CalculateHealthyPressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector);

/**
 * Function to calculate damaged pressures in the network
 * @param network
 * @param pressureSolutionVector
 */
void CalculateDamagedPressures(NetworkDescription& network, std::vector<double>& pressureSolutionVector);

/**
 * Function that initializes the rhs vector of the matrix equations.
 *
 * @param[in] network
 * @param[in,out] rhsSolutionVector
 */
void CreatePressureRhsVector(const NetworkDescription& network, std::vector<double>& rhsSolutionVector);