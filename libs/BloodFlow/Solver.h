#pragma once

#include <vector>

#include "Triplet.h"
#include "NetworkDescription.h"

void SolvePressures(const NetworkDescription& network, std::vector<Triplet>& tripletList, std::vector<double>& rhsOutputVector);