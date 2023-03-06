#pragma once

#include <set>
#include <string>
#include <vector>

#include <hdf5.h>

#include "EdgeFinderIO.h"
#include "EdgeFinderVessel.h"
#include "EdgeNode.h"
#include "PartitionInformation.h"


/**
 * Function that finds and labels the edge junctions in the network
 *
 * This function imports connected vessel data, not using the old data available
 * @param fileId
 * @param partitionVessels
 * @param partInfo
 * @param localEdgeNodes
 */

void FindAndLabelEdgeNodes(hid_t fileId,
    std::vector<EdgeFinderVessel>& partitionVessels,
    PartitionInformation& partInfo, std::set<EdgeNode>& localEdgeNodes);