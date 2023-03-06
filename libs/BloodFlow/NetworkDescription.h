#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "BloodFlowEdgeNode.h"
#include "BloodFlowVessel.h"

/**
 * Class to contain all items related to the vascular geomerty and partitioning
 *
 * I would probably rewrite this class to use more composition as it contain many many functions.
 */
class NetworkDescription
{
public:
    NetworkDescription();
    NetworkDescription(short localPartitionIndex);
    ~NetworkDescription();
    enum BoundrayConditionType
    {
        Pressure,
        Flow
    };

    friend std::ostream& operator<<(std::ostream& os, const NetworkDescription& rhs);

    long long getTotalNumberOfVessels() const
    {
        return mTotalNumberOfVessels;
    }
    void setTotalNumberOfVessels(long long val)
    {
        mTotalNumberOfVessels = val;
    }
    long long getLocalNumberOfVessels() const
    {
        return mNumberOFVesselsPerParition[mLocalPartitionIndex];
    }
    short getNumberOfPartitions() const
    {
        return mNumberOFVesselsPerParition.size();
    }
    short getLocalPartitionIndex() const
    {
        return mLocalPartitionIndex;
    }
    void setLocalPartitionIndex(short val)
    {
        mLocalPartitionIndex = val;
    }

    void setNumberOfVesselsInPartition(short partitionIndex, long long numVessels);
    long long getNumberOfVesselsInAPartition(short partitionIndex) const
    {
        return mNumberOFVesselsPerParition[partitionIndex];
    }
    std::vector<long long>& getNumberOfVesselsPerParition()
    {
        return mNumberOFVesselsPerParition;
    }
    const std::vector<long long>& getNumberOfVesselsPerParition() const
    {
        return mNumberOFVesselsPerParition;
    }

    void addRowBoundPairToVector(long long start, long long end);
    void addRowBoundPairToVector(std::pair<long long, long long>);
    void setRowBoundPairToVector(long long start, long long end, short partitionIndex);
    void setRowBoundPairToVector(std::pair<long long, long long>, short partitionIndex);
    std::pair<long long, long long>& getRowBoundPair(short partitionIndex)
    {
        return mRankRowBounds[partitionIndex];
    }
    const std::pair<long long, long long>& getRowBoundPair(
        short partitionIndex) const
    {
        return mRankRowBounds[partitionIndex];
    }
    long long getRankStartRow(short partitionIndex) const
    {
        return mRankRowBounds[partitionIndex].first;
    }
    long long getRankEndRow(short partitionIndex) const
    {
        return mRankRowBounds[partitionIndex].second;
    }
    long long getLocalNumRows() const
    {
        return getLocalRankEndRow() - getLocalRankStartRow() + 1;
    }
    long long getMatrixDimensions() const;

    std::pair<long long, long long>& getLocalRowBoundPair()
    {
        return mRankRowBounds[mLocalPartitionIndex];
    }
    const std::pair<long long, long long>& getLocalRowBoundPair() const
    {
        return mRankRowBounds[mLocalPartitionIndex];
    }
    long long getLocalRankStartRow() const
    {
        return mRankRowBounds[mLocalPartitionIndex].first;
    }
    long long getLocalRankEndRow() const
    {
        return mRankRowBounds[mLocalPartitionIndex].second;
    }
    long long determinePartitionForRow(long long rowIndex) const;

    void addEdgeNode(BloodFlowEdgeNode& edgeNode);
    const BloodFlowEdgeNode& getEdgeNode(long long nodeIndex) const
    {
        return *(mEdgeNodes.find(nodeIndex));
    }
    const std::set<BloodFlowEdgeNode>& getEdgeNodeSet() const
    {
        return mEdgeNodes;
    }
    int getEdgeNodeHostPartition(long long nodeIndex) const;
    int getNumLevels() const
    {
        return numLevels;
    }
    void setNumLevels(int val)
    {
        numLevels = val;
    }
    bool isEdgeNode(long long index) const;

    void addLocalVessel(BloodFlowVessel vessel);
    BloodFlowVessel& getLocalVessel(long long localIndex)
    {
        return mLocalVessels[localIndex];
    }
    const BloodFlowVessel& getLocalVessel(long long localIndex) const
    {
        return mLocalVessels[localIndex];
    }
    std::vector<BloodFlowVessel>& getLocalVesselVector()
    {
        return mLocalVessels;
    }
    const std::vector<BloodFlowVessel>& getLocalVesselVector() const
    {
        return mLocalVessels;
    }
    void clearLocalVesselVector()
    {
        mLocalVessels.clear();
    }

    long long getGlobalVesselStartIndex() const
    {
        return mGlobalVesselStartIndex;
    }
    void setGlobalVesselStartIndex(long long val)
    {
        mGlobalVesselStartIndex = val;
    }

    void addBoundaryConditionNode(long long nodeIndex, NetworkDescription::BoundrayConditionType type, double value);
    const std::pair<BoundrayConditionType, double>&
    getBoundaryConditionNodeInformation(long long index) const;
    std::pair<NetworkDescription::BoundrayConditionType, double>&
    getBoundaryConditionNodeInformation(long long index);
    std::map<long long, std::pair<BoundrayConditionType, double>>&
    getBoundaryConditionMap()
    {
        return mBoundaryConditionNodes;
    }
    const std::map<long long, std::pair<BoundrayConditionType, double>>&
    getBoundaryConditionMap() const
    {
        return mBoundaryConditionNodes;
    }
    bool isBoundaryConditionNode(long long index) const;
    const std::vector<std::pair<long long, long long>>& getRankBoundArray()
        const
    {
        return mRankRowBounds;
    }

    const std::vector<long long>& GetPermutedBcNodeNumbers() const
    {
        return mPermutedBcNodeNumbers;
    }
    std::vector<long long>& GetPermutedBcNodeNumbers()
    {
        return mPermutedBcNodeNumbers;
    }
    void AddPermutedBcIndex(long long index)
    {
        mPermutedBcNodeNumbers.push_back(index);
    }
    void ClearPermutedBcIndexes()
    {
        mPermutedBcNodeNumbers.clear();
    }

private:
    long long mTotalNumberOfVessels = 0;
    short mLocalPartitionIndex = 0;
    long long mGlobalVesselStartIndex = 0;
    int numLevels = 0;

    std::vector<long long> mNumberOFVesselsPerParition;
    std::vector<std::pair<long long, long long>> mRankRowBounds;
    std::set<BloodFlowEdgeNode> mEdgeNodes;
    std::map<long long, std::pair<BoundrayConditionType, double>> mBoundaryConditionNodes;
    std::vector<long long> mPermutedBcNodeNumbers;
    std::vector<BloodFlowVessel> mLocalVessels;
};