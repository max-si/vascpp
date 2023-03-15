#pragma once

#include <string>

#include "NetworkDescription.h"

void BloodFlowExe(std::string filename);

bool CheckHealthyFlowsAreCorrect(const NetworkDescription& network);