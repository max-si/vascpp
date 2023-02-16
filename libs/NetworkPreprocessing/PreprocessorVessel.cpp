
#include <algorithm>
#include <iostream>
#include <iterator>

#include "core/Timsort.h"
#include "PreprocessorVessel.h"

PreprocessorVessel::PreprocessorVessel()
  : Vessel()
  , partitionNumber(0)
{
}

// PreprocessorVessel::~PreprocessorVessel() {}

void PreprocessorVessel::setConnectedVessels(std::vector<long long>& input)
{
    connectedVessels.assign(input.begin(), input.end());
}

void PreprocessorVessel::assignConnectedVessels(
    std::vector<long long>::iterator& start,
    std::vector<long long>::iterator& end)
{
    connectedVessels.assign(start, end);
    gfx::timsort(connectedVessels.begin(), connectedVessels.end());
}

bool PreprocessorVessel::operator==(PreprocessorVessel& rhs)
{
    bool equal = true;
    equal = equal && (rhs.getPartitionNumber() == partitionNumber);
    equal = equal &&
        std::equal(connectedVessels.begin(), connectedVessels.end(),
            rhs.getConnectedVessels().begin());
    equal = equal && (rhs.getInitialArrayPosition() == initialArrayPosition);
    equal = equal && Vessel::operator==(rhs);

    return equal;
}
bool PreprocessorVessel::operator!=(PreprocessorVessel& rhs)
{
    return !(*this == rhs);
}

bool CompareVesselPartition(const PreprocessorVessel& lhs, const PreprocessorVessel& rhs)
{
    return lhs.getPartitionNumber() < rhs.getPartitionNumber();
}

bool CompareVesselNodesPartition(const PreprocessorVessel& lhs, const PreprocessorVessel& rhs)
{
    bool val;
    if (lhs.getPartitionNumber() < rhs.getPartitionNumber())
        return true;
    if (lhs.getPartitionNumber() > rhs.getPartitionNumber())
        return false;
    if (lhs.getPartitionNumber() == rhs.getPartitionNumber())
    {
        if (lhs.getNode1() < rhs.getNode1())
            return true;
        if (lhs.getNode1() > rhs.getNode1())
            return false;
        if (lhs.getNode2() < rhs.getNode2())
            return true;
        return false;
    }
    return false;
}

std::ostream& operator<<(std::ostream& os, const PreprocessorVessel& vessel)
{
    //os << vessel.getStartPoint() << "\t" << vessel.getEndPoint() << "\t" << vessel.getOrigRadius() << "\t" << vessel.getNode1() << "\t" << vessel.getNode2();
    os << "\n\tConnected Vessels (" << vessel.getNumConnectedVessels() << "): ";
    for (auto& x : vessel.getConnectedVessels())
        os << x << " ";
    os << "\n\tPartition Number: " << vessel.getPartitionNumber()
       << " \tOriginalIndex: " << vessel.getInitialArrayPosition();
    return os;
}