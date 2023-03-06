#pragma once

#include <vector>

#include "Triplet.h"
#include "NetworkDescription.h"

/**
 * Function to create a triplet list with healthy vessel conductances
 * @param[in] network
 * @param[out] tripletList
 */
void HealthyPressureListCreator(const NetworkDescription& network, std::vector<Triplet>& tripletList);

/**
 * Function to create a triplet list with damaged vessel conductances
 * @param[in] network
 * @param[out] tripletList
 */
void DamagedPressureListCreator(const NetworkDescription& network, std::vector<Triplet>& tripletList);