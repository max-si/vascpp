#pragma once

#include <string>
#include <vector>

#include <hdf5.h>

#include "BloodFlowVessel.h"
#include "NetworkDescription.h"

hid_t BloodFlowNetworkImporter(std::string fileName, NetworkDescription& networkDescription);