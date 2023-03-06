#include <cstdlib>

#include "PartitionInformation.h"

PartitionInformation::PartitionInformation() {}

PartitionInformation::~PartitionInformation()
{
    if (numConnectedVesselsArray)
    {
        delete[] numConnectedVesselsArray;
    }
    delete[] numVesselsPerPartition;
}

std::ostream& operator<<(std::ostream& os, PartitionInformation& partInfo)
{
    os << partInfo.numPartitions << " " << partInfo.numVesselsInPartition << " "
       << partInfo.startVesselIndex << " ";
    os << partInfo.numConnectedVessels << " "
       << partInfo.connectedVesselStartIndex << "\n\t";
    if (partInfo.numConnectedVesselsArray)
    {
        for (int i = 0; i < partInfo.numVesselsInPartition; i++)
            os << partInfo.numConnectedVesselsArray[i] << " ";
        os << "\n\t";
    }
    if (partInfo.numVesselsPerPartition)
    {
        for (int i = 0; i < partInfo.numPartitions + 1; i++)
            os << partInfo.numVesselsPerPartition[i] << " ";
    }
    return os;
}