#pragma once

#include <vector>

#include "NetworkDescription.h"

void HealthyFlowCalculations(NetworkDescription& network, std::vector<double>& pressureSolutionVector);

/**
 * \breif An support Function that adds Calculated healthy flow rates to each vessel in a vector
 *
 * A support function that adds the calculated blood flow rate to the vessel objects.
 * It is based on a simple for loop function
 *
 * @param[in,out] vessels A vector of blood vessel objects
 * @param[in] flowValues A vector of flow rates in the same sorted order as the vessel vector
 */
void AddHealthyFlowToVessels(
    std::vector<BloodFlowVessel>& vessels, std::vector<double>& flowValues);

/**
 * \breif An support Function that adds calculated damaged flow rates to each vessel in a vector
 *
 * A support function that adds the calculated blood flow rate to the vessel objects.
 * It is based on a simple for loop function
 *
 * @param[in,out] vessels A vector of blood vessel objects
 * @param[in] flowValues A vector of flow rates in the same sorted order as the vessel vector
 */
void AddDamagedFlowToVessels(
    std::vector<BloodFlowVessel>& vessels, std::vector<double>& flowValues);