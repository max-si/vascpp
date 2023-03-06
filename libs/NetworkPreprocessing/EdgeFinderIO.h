#pragma once

#include <set>
#include <vector>

#include <hdf5.h>

#include "EdgeFinderVessel.h"
#include "EdgeNode.h"
#include "PartitionInformation.h"

/**
 * Imports all the data necessary for the junction finding and constructs EdgeFinderVessels
 * @param fileId
 * @param partInfo
 * @param[out] vessels
 */
void ImportVesselsEdgeFinder(hid_t fileId, PartitionInformation& partInfo, std::vector<EdgeFinderVessel>& vessels);
/**
 * Imports node information from file and uses existing connected vessel data to create the vessel objects
 * @param fileId
 * @param partInfo
 * @param[out] vessels
 * @param connectedVesselData
 */
void ImportVesselsEdgeFinder(hid_t fileId, PartitionInformation& partInfo,
    std::vector<EdgeFinderVessel>& vessels,
    std::vector<long long> connectedVesselData);

void ExportEdgeNodesToFile(hid_t fileId, std::set<EdgeNode>& localEdgeNodes, PartitionInformation& partInfo);