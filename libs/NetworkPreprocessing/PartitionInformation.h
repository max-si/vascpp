#pragma once

#include <iostream>

/**
 * \breif Class that describes what a partition looklike in the network
 *
 * Stores all vessel data pertaining to vessel connections and number of vessels.
 */
class PartitionInformation
{
public:
    PartitionInformation();
    ~PartitionInformation();
    int partitionNum;
    int numPartitions;
    long long numVesselsInPartition;
    long long startVesselIndex;
    long long numConnectedVessels;
    long long connectedVesselStartIndex;
    short* numConnectedVesselsArray = nullptr;
    long long* numVesselsPerPartition = nullptr;

private:
    friend std::ostream& operator<<(
        std::ostream& os, PartitionInformation& partInfo);
};