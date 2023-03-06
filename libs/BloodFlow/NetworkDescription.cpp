#include <iostream>

#include "NetworkDescription.h"

NetworkDescription::NetworkDescription() {}

NetworkDescription::NetworkDescription(short localPartitionIndex)
  : mLocalPartitionIndex(localPartitionIndex) {}

NetworkDescription::~NetworkDescription() {}

void NetworkDescription::setNumberOfVesselsInPartition(
    short partitionIndex, long long numVessels)
{
    if (partitionIndex >= mNumberOFVesselsPerParition.size())
    {
        if (partitionIndex == 0)
        {
            mNumberOFVesselsPerParition.resize(1);
        }
        else
        {
            mNumberOFVesselsPerParition.resize(partitionIndex);
        }
    }
    mNumberOFVesselsPerParition[partitionIndex] = numVessels;
}

void NetworkDescription::addRowBoundPairToVector(long long start, long long end)
{
    mRankRowBounds.emplace_back(start, end);
}

void NetworkDescription::addRowBoundPairToVector(
    std::pair<long long, long long> val)
{
    mRankRowBounds.push_back(val);
}

void NetworkDescription::setRowBoundPairToVector(
    long long start, long long end, short partitionIndex)
{
    if (partitionIndex >= mRankRowBounds.size())
        mRankRowBounds.resize(partitionIndex);
    mRankRowBounds[partitionIndex].first = start;
    mRankRowBounds[partitionIndex].second = end;
}

void NetworkDescription::setRowBoundPairToVector(
    std::pair<long long, long long> val, short partitionIndex)
{
    if (partitionIndex >= mRankRowBounds.size())
        mRankRowBounds.resize(partitionIndex);
    mRankRowBounds[partitionIndex] = val;
}

long long NetworkDescription::getMatrixDimensions() const
{
    typedef std::pair<long long, long long> PairType;
    auto resultPair = *std::max_element(mRankRowBounds.begin(),
        mRankRowBounds.end(), [](const PairType& v1, const PairType& v2) {
            return v1.second < v2.second;
        });
    return resultPair.second + 1;
}

long long NetworkDescription::determinePartitionForRow(long long rowIndex) const
{
    for (int i = 0; i < mRankRowBounds.size(); ++i)
    {
        long long lowerBound = mRankRowBounds[i].first;
        long long upperBound = mRankRowBounds[i].second;
        if (rowIndex >= lowerBound && rowIndex <= upperBound)
            return i;
    }
    std::cout << "Partition Not Found for row index " << rowIndex << "\n";
    throw;
}

void NetworkDescription::addEdgeNode(BloodFlowEdgeNode& edgeNode)
{
    mEdgeNodes.insert(edgeNode);
}

int NetworkDescription::getEdgeNodeHostPartition(long long nodeIndex) const
{
    return getEdgeNode(nodeIndex).getHostPartition();
}

bool NetworkDescription::isEdgeNode(long long index) const
{
    return (mEdgeNodes.count(index));
}

void NetworkDescription::addLocalVessel(BloodFlowVessel vessel)
{
    mLocalVessels.push_back(vessel);
}

void NetworkDescription::addBoundaryConditionNode(
    long long nodeIndex, BoundrayConditionType type, double value)
{
    std::pair<BoundrayConditionType, double> myPair(std::make_pair(type, value));
    mBoundaryConditionNodes.insert(std::make_pair(nodeIndex, myPair));
}

std::pair<NetworkDescription::BoundrayConditionType, double>&
NetworkDescription::getBoundaryConditionNodeInformation(long long index)
{
    return mBoundaryConditionNodes[index];
}

const std::pair<NetworkDescription::BoundrayConditionType, double>&
NetworkDescription::getBoundaryConditionNodeInformation(long long index) const
{
    auto node = mBoundaryConditionNodes.find(index);
    return node->second;
}

bool NetworkDescription::isBoundaryConditionNode(long long index) const
{
    return (mBoundaryConditionNodes.count(index));
}

std::ostream& operator<<(std::ostream& os, const NetworkDescription& rhs)
{
    os << "Network Description Contents:\n"
       << "Total Number of Vessels: " << rhs.mTotalNumberOfVessels << '\n'
       << "Partition Index: " << rhs.mLocalPartitionIndex << '\n'
       << "Import Start Index: " << rhs.mGlobalVesselStartIndex << '\n'
       << "Number of Levels in Network: " << rhs.numLevels << '\n'
       << "Number of Vessels in Each Partition: ";
    const auto& numOfVesselsPerpartition = rhs.getNumberOfVesselsPerParition();
    for (int i = 0; i < numOfVesselsPerpartition.size(); i++)
        os << numOfVesselsPerpartition[i] << " ";
    os << '\n';
    if (rhs.mRankRowBounds.size())
    {
        os << "Row Bounds by Rank:\n";

        for (int i = 0; i < rhs.getNumberOfPartitions(); i++)
        {
            const auto& data = rhs.getRowBoundPair(i);
            os << "\t" << data.first << " " << data.second << '\n';
        }
    }

    if (rhs.mEdgeNodes.size())
    {
        const auto& thing = rhs.getEdgeNodeSet();
        os << "Edge Nodes:\n";
        for (auto itr = thing.begin(); itr != thing.end(); ++itr)
        {
            auto& datacopy = *itr;
            os << "\t" << datacopy << '\n';
        }
    }

    if (rhs.mBoundaryConditionNodes.size())
    {
        os << "Boundary Condition Nodes:\n";
        const auto& data2 = rhs.getBoundaryConditionMap();
        for (auto itr = data2.begin(); itr != data2.end(); ++itr)
        {
            os << "\t" << (*itr).first << ": " << (*itr).second.first << " - "
               << (*itr).second.second << '\n';
        }
        os << '\n';
    }

    os << "Permuted BC Node Indexes: " << rhs.mPermutedBcNodeNumbers.size()
       << " - ";
    if (rhs.mPermutedBcNodeNumbers.size())
    {
        for (long long i = 0; i < rhs.mPermutedBcNodeNumbers.size(); ++i)
        {
            os << "\t" << rhs.mPermutedBcNodeNumbers[i] << ' ';
        }
        os << '\n';
    }

    os << "Number Of Vessels: " << rhs.mLocalVessels.size() << '\n';
    if (rhs.mLocalVessels.size() < 70)
    {
        for (long long i = 0; i < rhs.mLocalVessels.size(); ++i)
        {
            os << "\t" << rhs.mLocalVessels[i] << '\n';
        }
    }

    return os;
}